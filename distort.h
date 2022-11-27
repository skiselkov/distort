/*
 * Copyright (c) 2022, Saso Kiselkov
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 * 1) Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2) Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 * 3) Neither the name of the copyright holder nor the names of its
 *    contributors may be used to endorse or promote products derived from
 *    this software without specific prior written permission.
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
 * IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
 * THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef	__DISTORT_H__
#define	__DISTORT_H__

#include <stdint.h>

#ifdef  __cplusplus
extern "C" {
#endif

typedef struct distort_s distort_t;

/*
 * Initializer and finalizer functions. Call to create & destroy context.
 * Sample rate must be one of: 44100 or 48000.
 */
distort_t *distort_init(unsigned sample_rate);
void distort_fini(distort_t *dis);

/*
 * Applies distortion to an input buffer. The buffer is replaced by the
 * distorted audio samples. If there isn't enough data in the input buffer
 * or cached data in the `dis' state object, the buffer will have silence
 * prepended to it to avoid outputting undistorted audio.
 *
 * Argument description:
 *
 * `dis': distortion control object previously created using distort_init().
 * `samples': an array of 16-bit signed PCM samples (single channel only!).
 * `num_samples': number of samples contained in `samples'.
 * `amplify': voice signal amplification level. You can pass 1.0 in here
 *	to not have the signal amplified at all, or a value less than 1.0
 *	to suppress the signal. This is not applied to the background noise
 *	generator, so you simulate a fading volume on the signal portion
 *	while keeping the noise volume constant.
 * `noise_level': background noise level from 0.0 to 1.0 (nothing but noise
 *	in the output). Typical values would be 0.02 for little background
 *	noise, 0.2 for moderate background noise, to 0.6 for heavy background
 *	noise (makes transmission almost unreadable).
 */
void distort(distort_t *dis, int16_t *samples, size_t num_samples,
    double amplify, double noise_level);

/*
 * Clears any cached buffers from the distortion object. Call this between
 * audio transmissions to avoid any partially completed chunks being played
 * at the start of an unrelated transmission.
 */
void distort_clear_buffers(distort_t *dis);

#ifdef  __cplusplus
}
#endif

#endif	/* __DISTORT_H__ */
