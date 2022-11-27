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
			buf = (uint8_t *)realloc(buf, bufcap);
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
