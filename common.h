/*
 * Copyright (C) 2013-2015 Intel Corporation
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 *
 */

#define TEMP_RECORD_FILE_NAME "/tmp/test.wav"

#define COMPOSE(a, b, c, d) ((a) | ((b)<<8) | ((c)<<16) | ((d)<<24))
#define WAV_RIFF			COMPOSE('R', 'I', 'F', 'F')
#define WAV_WAVE			COMPOSE('W', 'A', 'V', 'E')
#define WAV_FMT				COMPOSE('f', 'm', 't', ' ')
#define WAV_DATA			COMPOSE('d', 'a', 't', 'a')
#define WAV_FORMAT_PCM			1	/* PCM WAVE file encoding */

#define MAX_NUMBER_OF_CHANNELS		2
#define MAX_NB_OF_PEAK			10
#define MAX_NB_OF_FRAMES		(10 * 1024 * 1024)
#define PLAYBACK_TIME_BEFORE_CAPTURE	500	/* Given in ms */

#define ENOPEAK				100
#define EONLYDC				101
#define EBADPEAK			102

#define DC_THRESHOLD			7.01

#define FOUND_DC			(1<<1)
#define FOUND_WRONG_PEAK		(1<<0)

#define CHANNEL_MAX			2
#define CHANNEL_MIN			1

/* memory error message */
#define E_MALLOC	"error: Allocate memory buffer: "
/* file error message */
#define E_OPENFILEP	"error: Open file for playback: "
#define E_OPENFILEC	"error: Open file for capture: "
#define E_READFILE	"error: Read file: "
#define E_WRITEFILE	"error: Write file: "
#define E_SEEKFILE	"error: Seek file: "
#define E_FILECONTENT	"error: File content: "
/* thread error message */
#define E_NEWTHREADP	"error: Can't create playback thread: "
#define E_NEWTHREADC	"error: Can't create capture thread: "
#define E_JOINTHREADP	"error: Can't join playback thread: "
#define E_JOINTHREADC	"error: Can't join capture thread: "
#define E_EXITTHREADP	"error: Exit playback thread fail: "
#define E_EXITTHREADC	"error: Exit capture thread fail: "
/* pcm device error message */
#define E_OPENPCMP	"error: Open PCM playback device: "
#define E_OPENPCMC	"error: Open PCM capture device: "
#define E_READPCM	"error: Read PCM device: "
#define E_WRITEPCM	"error: Write PCM device: "
#define E_NOPCMP	"error: No PCM device for playback: "
#define E_NOPCMC	"error: No PCM device for capture: "
#define E_SETDEV	"error: Set parameter to device: "
#define E_GETDEV	"error: Get parameter from device: "
/* parameter error message */
#define E_PARAMS	"error: Invalid parameter: "
/* error type string */
#define S_DEFAULT	"default params: "
#define S_HWPARAMS	"hw params: "
#define S_CHANNELS	"channel number: "
#define S_PCMFORMAT	"PCM format: "
#define S_ACCESS	"access type: "
#define S_SAMPLERATE	"sample rate: "
#define S_BUFFERTIME	"buffer time: "
#define S_BUFFERSIZE	"buffer size: "
#define S_PERIODTIME	"period time: "
#define S_PERIODSIZE	"period size: "
#define S_UNDERRUN	"underrun: "
#define S_OVERRUN	"overrun: "

#define loge(key, ...)	do {		\
	fprintf(stderr, key);		\
	fprintf(stderr, ##__VA_ARGS__);	\
	fprintf(stderr, "\n");		\
} while (0)

struct wav_header {
	unsigned int magic; /* 'RIFF' */
	unsigned int length; /* file len */
	unsigned int type; /* 'WAVE' */
};

struct wav_chunk_header {
	unsigned int type; /* 'data' */
	unsigned int length; /* sample count */
};

struct wav_fmt {
	unsigned int magic; /* 'FMT '*/
	unsigned int fmt_size; /* 16 or 18 */
	unsigned short format; /* see WAV_FMT_* */
	unsigned short channels;
	unsigned int sample_rate; /* Frequency of sample */
	unsigned int bytes_p_second;
	unsigned short blocks_align; /* sample size; 1 or 2 bytes */
	unsigned short sample_length; /* 8, 12 or 16 bit */
};

struct chunk_fmt {
	unsigned short format; /* see WAV_FMT_* */
	unsigned short channels;
	unsigned int sample_rate; /* Frequency of sample */
	unsigned int bytes_p_second;
	unsigned short blocks_align; /* sample size; 1 or 2 bytes */
	unsigned short sample_length; /* 8, 12 or 16 bit */
};

struct wav_container {
	struct wav_header header;
	struct wav_fmt format;
	struct wav_chunk_header chunk;
};

struct bat;

struct play_cap {
	char *device;
	unsigned int card_tiny;
	unsigned int device_tiny;
	char *file;
	bool single;
	void *(*fct)(struct bat *);
};

struct bat {
	unsigned int rate;	/* sampling rate */
	int channels;		/* nb of channels */
	int frames;			/* nb of frames */
	int frame_size;		/* size of frame */
	int sample_size;	/* size of sample */

	float sigma_k;		/* threshold for peak detection */
	float target_freq[MAX_NUMBER_OF_CHANNELS];

	int sinus_duration;	/* */
	char *narg;			/* */

	struct play_cap playback;
	struct play_cap capture;

	unsigned int periods_played;
	unsigned int periods_total;
	bool period_limit;
	FILE *fp;

	double (*convert_sample_to_double)(void *, int);
	void (*convert_float_to_sample)(float, void *);

	void *buf;			/* PCM Buffer */

	bool local;			/* true for internal test */
	bool tinyalsa;		/* true if tinyalsa is used */
};

struct analyze {
	void *buf;
	double *in;
	double *out;
	double *mag;
};

void close_file(FILE *);
void destroy_mem(void *);

void prepare_wav_info(struct wav_container *, struct bat *);
int read_wav_header(struct bat *, char *, FILE *, bool);
void generate_sine_wave(struct bat *, int , void *, int);
