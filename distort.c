
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "distort.h"
#include "kiss_fft.h"

/*
 * Chunk length expressed as a fraction of a second (1 / TIME_QUANTUM).
 * This must evenly divide the sample rate, otherwise bad things happen.
 */
#define	TIME_QUANTUM		30
/*
 * Compressor target signal level. Signals weaker than this get amplified.
 * Signals stronger than this get suppressed.
 */
#define	COMPR_TGT		0.66	/* dimensionless */
/*
 * Compressor signal minimum energy level. This limits by how much we can
 * amplify sounds in th compressor (0.6 = max amplification 40%).
 */
#define	COMPR_MIN_ENERGY	0.5	/* dimensionless */

/*
 * How quickly the noise level randomizer jumps around as fraction of
 * time quantum.
 */
#define	NOISE_RAND_RATE		2.0	/* FILTER_IN rate arg */

struct distort_s {
	unsigned	srate;
	double		compr_energy;
	double		noise_level;
	double		noise_level_cur;
	double		amplify;

	int		chunk_a_b;

	int16_t		*inbuf;
	size_t		inbuf_fill;
	size_t		inbuf_cap;

	int16_t		*outbuf;
	size_t		outbuf_fill;
	size_t		outbuf_fill_act;
	size_t		outbuf_cap;
	/* Used by distortion algorithm */
	int16_t		*tmpbuf;

	/* Used by EQ */
	kiss_fft_cpx	*fin;
	kiss_fft_cpx	*fout;

	kiss_fft_cfg	cfg;
	kiss_fft_cfg	cfg_inv;
};

static void distort_EQ(distort_t *dis, const int16_t *samples,
    int16_t *out_samples, size_t num_samples);
#ifdef	XXX
static void distort_compressor(distort_t *dis, int16_t *samples,
    size_t num_samples);
#endif	/* XXX */
static void distort_noise_randomize(distort_t *dis, size_t num_samples);
static void distort_process(distort_t *dis);
static void distort_impl(distort_t *dis, const int16_t *in_samples,
    int16_t *out_samples, size_t num_samples, double amplify,
    double noise_level);

/*
 * A generic CRC64 implementation from OpenSolaris. Be sure to call
 * crc64_init before using it. Then just call crc64() and pass it the
 * data you want checksummed. Also includes a fast & light-weight
 * portable pseudo random number generator.
 */

#define CRC64_POLY      0xC96C5795D7870F42ULL   /* ECMA-182, reflected form */

static uint64_t crc64_table[256];
static uint64_t rand_seed = 0;

static void
crc64_init(void)
{
	int i;

	for (i = 0; i < 256; i++) {
		uint64_t *ct;
		int j;
		for (ct = crc64_table + i, *ct = i, j = 8; j > 0; j--)
			*ct = (*ct >> 1) ^ (-(*ct & 1) & CRC64_POLY);
	}
}

/*
 * Computes the CRC64 checksum of a block of input data.
 */
static uint64_t
crc64(const void *input, size_t sz)
{
	const uint8_t *in_bytes = input;
	uint64_t crc = -1ULL;
	size_t i;

	assert(crc64_table[128] == CRC64_POLY);
	for (i = 0; i < sz; i++)
		crc = (crc >> 8) ^ crc64_table[(crc ^ in_bytes[i]) & 0xFF];

	return (crc);
}

/*
 * Initializes the CRC64-based pseudo random number generator. Pass in some
 * random seed (e.g. current microclock() usually does nicely). Obviously
 * you only want to call this once in your app.
 */
static void
crc64_srand(uint64_t seed)
{
	rand_seed = seed;
}

/*
 * Grabs a random 64-bit number from the PRNG. This function isn't
 * thread-safe, so take care not to rely on its output being super-duper
 * unpredictable in multi-threaded apps.
 */
static uint64_t
crc64_rand(void)
{
	rand_seed = crc64(&rand_seed, sizeof (rand_seed));
	return (rand_seed);
}

static inline int
clampi(int x, int min_val, int max_val)
{
	if (x < min_val)
		return (min_val);
	if (x > max_val)
		return (max_val);
	return (x);
}

#define	NULL_VECT2		((vect2_t){NAN, NAN})
#define	IS_NULL_VECT(a)		(isnan((a).x))
#define	VECT2(x, y)		((vect2_t){(x), (y)})
typedef struct vect2_s {
	double	x;
	double	y;
} vect2_t;

/*
 * Interpolates a linear function defined by two points.
 *
 * @param x Point who's 'y' value we're looking for on the function.
 * @param x1 First reference point's x coordinate.
 * @param y1 First reference point's y coordinate.
 * @param x2 Second reference point's x coordinate.
 * @param y2 Second reference point's y coordinate.
 */
static double
fx_lin(double x, double x1, double y1, double x2, double y2)
{
        return (((x - x1) / (x2 - x1)) * (y2 - y1) + y1);
}

/*
 * Multi-segment version of fx_lin. The segments are defined as a series of
 * x-y coordinate points (list terminated with a NULL_VECT2). The list must
 * contain AT LEAST 2 points. The value of 'x' is then computed using the
 * fx_lin function from the appropriate segment. If 'x' falls outside of the
 * curve range, the `extrapolate' argument controls behavior. If extrapolate
 * is 1, the nearest segment is extrapolated to the value of 'x'.
 * Otherwise the function returns NAN.
 */
static double
fx_lin_multi(double x, const vect2_t *points, int extrapolate)
{
	assert(!IS_NULL_VECT(points[0]));
	assert(!IS_NULL_VECT(points[1]));

	for (;;) {
		vect2_t p1 = points[0], p2 = points[1];

		assert(p1.x < p2.x);

		if (x < p1.x) {
			/* X outside of range to the left */
			if (extrapolate)
				return (fx_lin(x, p1.x, p1.y, p2.x, p2.y));
			break;
		}
		/* X in range of current segment */
		if (x <= p2.x)
			return (fx_lin(x, p1.x, p1.y, p2.x, p2.y));
		/* X outside of range to the right */
		if (IS_NULL_VECT(points[2])) {
			if (extrapolate)
				return (fx_lin(x, p1.x, p1.y, p2.x, p2.y));
			break;
		}

		points++;
	}

	return (NAN);
}

/*
 * Weighted avg, 'w' is weight fraction from 0.0 = all of x to 1.0 = all of y.
 */
static inline double
wavg(double x, double y, double w)
{
	return (x + (y - x) * w);
}

/*
 * Provides a gradual method of integrating an old value until it approaches
 * a new target value. This is used in iterative processes by calling the
 * FILTER_IN macro repeatedly a certain time intervals (d_t = delta-time).
 * As time progresses, old_val will gradually be made to approach new_val.
 * The lag serves to make the approach slower or faster (e.g. a value of
 * '2' and d_t in seconds makes old_val approach new_val with a ramp that
 * is approximately 2 seconds long).
 */
#define FILTER_IN(old_val, new_val, d_t, lag) \
	do { \
		double o = (old_val); \
		double n = (new_val); \
		(old_val) += (n - o) * ((d_t) / (lag)); \
		/* Prevent an overshoot */ \
		if ((o < n && (old_val) > n) || \
		    (o > n && (old_val) < n)) \
			(old_val) = n; \
	} while (0)

distort_t *
distort_init(unsigned sample_rate)
{
	distort_t *dis = calloc(1, sizeof (*dis));
	size_t chunksz;

	/* Init shared state */
	crc64_init();
	crc64_srand(time(NULL));

	/* Only supported sample rates ATM */
	assert(sample_rate == 44100 || sample_rate == 48000);

	dis->srate = sample_rate;
	dis->compr_energy = 1.0;
	dis->amplify = 1.0;

	chunksz = dis->srate / TIME_QUANTUM;
	dis->tmpbuf = calloc(chunksz, sizeof (*dis->tmpbuf));
	dis->fin = calloc(chunksz, sizeof (*dis->fin));
	dis->fout = calloc(chunksz, sizeof (*dis->fout));
	dis->cfg = kiss_fft_alloc(chunksz, 0, NULL, NULL);
	dis->cfg_inv = kiss_fft_alloc(chunksz, 1, NULL, NULL);

	return (dis);
}

void
distort_fini(distort_t *dis)
{
	free(dis->inbuf);
	free(dis->outbuf);
	free(dis->tmpbuf);

	free(dis->fin);
	free(dis->fout);

	free(dis->cfg);
	free(dis->cfg_inv);

	free(dis);
}

void
distort(distort_t *dis, int16_t *samples, size_t num_samples,
    double amplify, double noise_level)
{
	int16_t *out_samples = calloc(num_samples, sizeof (*out_samples));

	distort_impl(dis, samples, out_samples, num_samples, amplify,
	    noise_level);
	memcpy(samples, out_samples, num_samples * sizeof (*samples));
	free(out_samples);
}

void
distort_clear_buffers(distort_t *dis)
{
	dis->inbuf_fill = 0;
	dis->outbuf_fill = 0;
	dis->outbuf_fill_act = 0;
}

static void
distort_impl(distort_t *dis, const int16_t *in_samples, int16_t *out_samples,
    size_t num_samples, double amplify, double noise_level)
{
	size_t i = 0;
	size_t chunksz = dis->srate / TIME_QUANTUM;
	size_t avail;

	dis->amplify = amplify;
	dis->noise_level = noise_level;

	/* Grow our input buffer if necessary */
	if (dis->inbuf_fill + num_samples > dis->inbuf_cap) {
		dis->inbuf_cap = dis->inbuf_fill + num_samples;
		dis->inbuf = realloc(dis->inbuf,
		    sizeof (*dis->inbuf) * dis->inbuf_cap);
	}
	/* Stick newly incoming data on the end of our input buffer */
	memcpy(&dis->inbuf[dis->inbuf_fill], in_samples,
	    sizeof (*in_samples) * num_samples);
	dis->inbuf_fill += num_samples;

	distort_process(dis);

	/*
	 * Populate output with whatever we have accumulated in our output
	 * buffer. If we don't have enough, pad the leading portion with
	 * silence.
	 */
	if (dis->outbuf_fill < num_samples + chunksz) {
		/*
		 * Still in the growing phase of the buffer, wait until there
		 * is enough data there to have a little extra.
		 */
		avail = (dis->outbuf_fill >= 2 * chunksz ?
		    dis->outbuf_fill - 2 * chunksz : 0);
	} else {
		/*
		 * Now with the incoming data, we should always have enough,
		 * so start handing it out.
		 */
		avail = (dis->outbuf_fill >= chunksz ?
		    dis->outbuf_fill - chunksz : 0);
	}
	if (num_samples > avail) {
		i = num_samples - avail;
		memset(out_samples, 0, sizeof (*out_samples) * i);
	}
	if (i < num_samples) {
		size_t wanted = num_samples - i;
		size_t to_copy = (wanted < avail ? wanted : avail);

		memcpy(&out_samples[i], dis->outbuf,
		    to_copy * sizeof (*dis->outbuf));

		/* shift remainder in output buffer forward */
		memmove(dis->outbuf, &dis->outbuf[to_copy],
		    sizeof (int16_t) * (dis->outbuf_fill_act - to_copy));
		dis->outbuf_fill -= to_copy;
		dis->outbuf_fill_act -= to_copy;
	}
}

static void
distort_process(distort_t *dis)
{
	size_t consumed;
	size_t chunksz = dis->srate / TIME_QUANTUM;

	for (consumed = 0; consumed + chunksz < dis->inbuf_fill;
	    consumed += chunksz / 2) {
		int16_t *samples = &dis->inbuf[consumed];

		/*
		 * Full chunk processing only happens on odd chunk.
		 * This is when we want to apply the compressor & noise
		 * level randomizer (don't want to have those differ
		 * between EQ passes).
		 */
		if (dis->chunk_a_b == 0) {
#ifdef	XXX
			distort_compressor(dis, samples, chunksz);
#endif
			distort_noise_randomize(dis, chunksz);
		}

		distort_EQ(dis, samples, dis->tmpbuf, chunksz);

		/* now mix in the latter half of the previous buffer */
		if (dis->outbuf_fill_act != 0) {
			size_t i;

			for (i = 0; i < chunksz / 2; i++) {
				int16_t oldval =
				    dis->outbuf[dis->outbuf_fill + i];
				int16_t newval = dis->tmpbuf[i];

				dis->tmpbuf[i] = wavg(oldval, newval,
				    i / (chunksz * 0.5));
			}
		}

		/*
		 * Grow output buffer if necessary.
		 */
		if (dis->outbuf_fill + chunksz > dis->outbuf_cap) {
			dis->outbuf_cap += chunksz;
			dis->outbuf = realloc(dis->outbuf,
			    sizeof (*dis->outbuf) * dis->outbuf_cap);
		}
		memcpy(&dis->outbuf[dis->outbuf_fill],
		    dis->tmpbuf, sizeof (*dis->tmpbuf) * chunksz);
		dis->outbuf_fill += chunksz / 2;
		dis->outbuf_fill_act = dis->outbuf_fill + chunksz / 2;

		dis->chunk_a_b = !dis->chunk_a_b;
	}

	/* Shift uncomsumed data in the input buffer forward. */
	if (consumed > 0) {
		assert(consumed < dis->inbuf_fill);
		memmove(dis->inbuf, &dis->inbuf[consumed],
		    sizeof (int16_t) * (dis->inbuf_fill - consumed));
		dis->inbuf_fill -= consumed;
	}
}

#ifdef	XXX
static void
distort_compressor(distort_t *dis, int16_t *samples, size_t num_samples)
{
	size_t i;

	for (i = 0; i < num_samples; i++) {
		double sample = (samples[i] >= 0 ? samples[i] : -samples[i]);
		double e = sample / (INT16_MAX * 1.0);

		e = (e > COMPR_MIN_ENERGY ? e : COMPR_MIN_ENERGY);
		if (e >= dis->compr_energy)
			FILTER_IN(dis->compr_energy, 1, 1.0, 100.0);
		else
			FILTER_IN(dis->compr_energy, 0, 1.0, 1000.0);
		samples[i] = clampi(sample / dis->compr_energy,
		    INT16_MIN, INT16_MAX);
	}
}
#endif	/* XXX */

static void
distort_noise_randomize(distort_t *dis, size_t num_samples)
{
	/*
	 * Noise level oscillates between 0.5x - 1.5x requested level.
	 */
	double noise_level_tgt = dis->noise_level *
	    (1 + ((crc64_rand() / (double)UINT64_MAX - 0.5)));

	FILTER_IN(dis->noise_level_cur, noise_level_tgt, num_samples,
	    num_samples * NOISE_RAND_RATE);
}

static void
distort_EQ(distort_t *dis, const int16_t *samples, int16_t *out_samples,
    size_t num_samples)
{
#define	HZ2SLOT(hz)	(((hz) / (double)dis->srate) * num_samples)
#define	CENTER_AMPLIFY	1.6
	const vect2_t clamp_curve[] = {
	    VECT2(HZ2SLOT(0), 0),
	    VECT2(HZ2SLOT(250), 0),
	    VECT2(HZ2SLOT(300), num_samples),
	    VECT2(HZ2SLOT(1700), CENTER_AMPLIFY * num_samples),
	    VECT2(HZ2SLOT(3000), num_samples),
	    VECT2(HZ2SLOT(3500), 0),
	    VECT2(HZ2SLOT(44500), 0),
	    VECT2(HZ2SLOT(45000), num_samples),
	    VECT2(HZ2SLOT(46300), CENTER_AMPLIFY * num_samples),
	    VECT2(HZ2SLOT(47700), num_samples),
	    VECT2(HZ2SLOT(47750), 0),
	    VECT2(HZ2SLOT(48000), 0),
	    NULL_VECT2
	};
	size_t i;

	for (i = 0; i < num_samples; i++) {
		/* apply noise to the input */
		int16_t rand_sample = (int16_t)crc64_rand();
		dis->fin[i].r = dis->amplify * samples[i] +
		    rand_sample * dis->noise_level_cur;
		dis->fin[i].i = 0;
	}
	/*
	 * Transform from time domain to frequency domain.
	 */
	kiss_fft(dis->cfg, dis->fin, dis->fout);
	for (i = 0; i < num_samples; i++) {
		double scale = fx_lin_multi(i, clamp_curve, 1);
		/*
		 * Rescale FFT factors to apply the EQ in frequency domain.
		 */
		dis->fout[i].r = dis->fout[i].r * scale;
		dis->fout[i].i = dis->fout[i].i * scale;
	}
	/*
	 * Inverse transform from frequency domain back to time domain.
	 */
	kiss_fft(dis->cfg_inv, dis->fout, dis->fin);

	for (i = 0; i < num_samples; i++) {
		/*
		 * Sometimes the values just hit the max/min values,
		 * so clamp it properly to avoid overflow/underflow.
		 */
		out_samples[i] = clampi(dis->fin[i].r, INT16_MIN, INT16_MAX);
	}
}
