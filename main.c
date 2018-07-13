#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include <acfutils/riff.h>

#include "distort.h"

#define	READSZ	32768
#define	WAVE_ID	FOURCC("WAVE")
#define	DATA_ID	FOURCC("data")

int
main(void)
{
	uint8_t *buf = NULL;
	size_t bufcap = 0, bufsz = 0;
	riff_chunk_t *riff, *data;
	distort_t *dis;
	size_t i;

	while (!feof(stdin)) {
		if (bufsz + READSZ > bufcap) {
			bufcap += READSZ;
			buf = realloc(buf, bufcap);
		}
		bufsz += fread(&buf[bufsz], 1, bufcap - bufsz, stdin);
		if (ferror(stdin)) {
			perror("Error reading stdin");
			return (1);
		}
	}

	riff = riff_parse(WAVE_ID, buf, bufsz);
	if (riff == NULL) {
		fprintf(stderr, "Malformed WAV: failed to parse RIFF\n");
		return (1);
	}
	data = riff_find_chunk(riff, DATA_ID, 0);
	if (data == NULL) {
		fprintf(stderr, "Malformed WAV: no data\n");
		return (1);
	}

	dis = distort_init(48000);
#define	CHUNKSZ	6783
	for (i = 0; i < data->datasz / 2; i += CHUNKSZ)
		distort(dis, &((int16_t *)data->data)[i], CHUNKSZ, 1.0, 0.0);
	fwrite(buf, 1, bufsz, stdout);
	distort_fini(dis);

	free(buf);

	return (0);
}
