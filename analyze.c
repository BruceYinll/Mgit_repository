#include <stdio.h> 
21 #include <stdlib.h> 
22 #include <errno.h> 
23 #include <stdbool.h> 
24 #include <stdint.h> 
25 
 
26 #include <math.h> 
27 #include <fftw3.h> 
28 
 
29 #include "common.h" 
30 
 
31 /** 
32  * Convert from sample size to double 
33  */ 
34 static int convert(struct bat *bat, struct analyze *a) 
35 { 
36 	void *s = a->buf; 
37 	int i; 
38 
 
39 	for (i = 0; i < bat->frames; i++) 
40 		a->in[i] = bat->convert_sample_to_double(s, i); 
41 
 
42 	return 0; 
43 } 
44 /** 
45  * 
46  * @return 0 if peak detected at right frequency, 1 if peak detected somewhere else 
47  *         2 if DC detected 
48  */ 
49 int check_peak(struct bat *bat, struct analyze *a, int end, int peak, float hz, 
50 		float mean, float p, int channel, int start) 
51 { 
52 	int ret; 
53 	float hz_peak = (float) (peak) * hz; 
54 	float delta_rate = 0.005 * bat->target_freq[channel]; 
55 	float delta_HZ = 1.0; 
56 	float tolerance = (delta_rate > delta_HZ) ? delta_rate : delta_HZ; 
57 
 
58 	printf("Detected peak at %2.2f Hz of %2.2f dB\n", hz_peak, 
59 			10.0 * log10(a->mag[peak] / mean)); 
60 	printf(" Total %3.1f dB from %2.2f to %2.2f Hz\n", 
61 			10.0 * log10(p / mean), start * hz, 
62 			end * hz); 
63 
 
64 	if (hz_peak < DC_THRESHOLD) { 
65 		fprintf(stdout, 
66 				" WARNING: Found low peak %2.2f Hz, very close to DC\n", 
67 				hz_peak); 
68 		ret = FOUND_DC; 
69 	} else if (hz_peak < bat->target_freq[channel] - tolerance) { 
70 		printf(" FAIL: Peak freq too low %2.2f Hz\n", hz_peak); 
71 		ret = FOUND_WRONG_PEAK; 
72 	} else if (hz_peak > bat->target_freq[channel] + tolerance) { 
73 		fprintf(stdout, 
74 				" FAIL: Peak freq too high %2.2f Hz\n", 
75 				hz_peak); 
76 		ret = FOUND_WRONG_PEAK; 
77 	} else { 
78 		printf(" PASS: Peak detected at target frequency\n"); 
79 		ret = 0; 
80 	} 
81 
 
82 	return ret; 
83 } 
84 
 
85 /** 
86  * Search for main frequencies in fft results and compare it to target 
87  */ 
88 static int check(struct bat *bat, struct analyze *a, int channel) 
89 { 
90 	float hz = 1.0 / ((float) bat->frames / (float) bat->rate); 
91 	float mean = 0.0, t, sigma = 0.0, p = 0.0; 
92 	int i, start = -1, end = -1, peak = 0, signals = 0; 
93 	int ret = 0, N = bat->frames / 2; 
94 
 
95 	/* calculate mean */ 
96 	for (i = 0; i < N; i++) 
97 		mean += a->mag[i]; 
98 	mean /= (float) N; 
99 
 
100 	/* calculate standard deviation */ 
101 	for (i = 0; i < N; i++) { 
102 		t = a->mag[i] - mean; 
103 		t *= t; 
104 		sigma += t; 
105 	} 
106 	sigma /= (float) N; 
107 	sigma = sqrtf(sigma); 
108 
 
109 	/* clip any data less than k sigma + mean */ 
110 	for (i = 0; i < N; i++) { 
111 		if (a->mag[i] > mean + bat->sigma_k * sigma) { 
112 
 
113 			/* find peak start points */ 
114 			if (start == -1) { 
115 				start = peak = end = i; 
116 				signals++; 
117 			} else { 
118 				if (a->mag[i] > a->mag[peak]) 
119 					peak = i; 
120 				end = i; 
121 			} 
122 			p += a->mag[i]; 
123 		} else if (start != -1) { 
124 			/* Check if peak is as expected */ 
125 			ret |= check_peak(bat, a, end, peak, hz, mean, 
126 					p, channel, start); 
127 			end = start = -1; 
128 			if (signals == MAX_NB_OF_PEAK) 
129 				break; 
130 		} 
131 	} 
132 	if (signals == 0) 
133 		ret = -ENOPEAK; /* No peak detected */ 
134 	else if ((ret == FOUND_DC) && (signals == 1)) 
135 		ret = -EONLYDC; /* Only DC detected */ 
136 	else if ((ret & FOUND_WRONG_PEAK) == FOUND_WRONG_PEAK) 
137 		ret = -EBADPEAK; /* Bad peak detected */ 
138 	else 
139 		ret = 0; /* Correct peak detected */ 
140 
 
141 	printf("Detected at least %d signal(s) in total\n", signals); 
142 
 
143 	return ret; 
144 } 
145 
 
146 static void calc_magnitude(struct bat *bat, struct analyze *a, int N) 
147 { 
148 	double r2, i2; 
149 	int i; 
150 
 
151 	for (i = 1; i < N / 2; i++) { 
152 		r2 = a->out[i] * a->out[i]; 
153 		i2 = a->out[N - i] * a->out[N - i]; 
154 
 
155 		a->mag[i] = sqrt(r2 + i2); 
156 	} 
157 	a->mag[0] = 0.0; 
158 } 
159 
 
160 static int find_and_check_harmonics(struct bat *bat, struct analyze *a, 
161 		int channel) 
162 { 
163 	fftw_plan p; 
164 	int ret = -ENOMEM, N = bat->frames; 
165 
 
166 	/* Allocate FFT buffers */ 
167 	a->in = (double *) fftw_malloc(sizeof(double) * bat->frames); 
168 	if (a->in == NULL) 
169 		goto out1; 
170 
 
171 	a->out = (double *) fftw_malloc(sizeof(double) * bat->frames); 
172 	if (a->out == NULL) 
173 		goto out2; 
174 
 
175 	a->mag = (double *) fftw_malloc(sizeof(double) * bat->frames); 
176 	if (a->mag == NULL) 
177 		goto out3; 
178 
 
179 	/* create FFT plan */ 
180 	p = fftw_plan_r2r_1d(N, a->in, a->out, FFTW_R2HC, 
181 			FFTW_MEASURE | FFTW_PRESERVE_INPUT); 
182 	if (p == NULL) 
183 		goto out4; 
184 
 
185 	/* convert source PCM to doubles */ 
186 	ret = convert(bat, a); 
187 	if (ret != 0) 
188 		goto out4; 
189 
 
190 	/* run FFT */ 
191 	fftw_execute(p); 
192 
 
193 	/* FFT out is real and imaginary numbers - calc magnitude for each */ 
194 	calc_magnitude(bat, a, N); 
195 
 
196 	/* check data */ 
197 	ret = check(bat, a, channel); 
198 
 
199 	fftw_destroy_plan(p); 
200 
 
201 out4: 
202 	fftw_free(a->mag); 
203 out3: 
204 	fftw_free(a->out); 
205 out2: 
206 	fftw_free(a->in); 
207 out1: 
208 	return ret; 
209 } 
210 
 
211 /** 
212  * Convert interleaved samples from channels in samples from a single channel 
213  */ 
214 static int reorder_data(struct bat *bat) 
215 { 
216 	char *p, *new_bat_buf; 
217 	int ch, i, j; 
218 
 
219 	if (bat->channels == 1) 
220 		return 0; /* No need for reordering */ 
221 
 
222 	p = malloc(bat->frames * bat->frame_size); 
223 	new_bat_buf = p; 
224 	if (p == NULL) 
225 		return -ENOMEM; 
226 
 
227 	for (ch = 0; ch < bat->channels; ch++) { 
228 		for (j = 0; j < bat->frames; j++) { 
229 			for (i = 0; i < bat->sample_size; i++) { 
230 				*p++ = ((char *) (bat->buf))[j * bat->frame_size 
231 						+ ch * bat->sample_size + i]; 
232 			} 
233 		} 
234 	} 
235 
 
236 	free(bat->buf); 
237 	bat->buf = new_bat_buf; 
238 
 
239 	return 0; 
240 } 
241 
 
242 int analyze_capture(struct bat *bat) 
243 { 
244 	int ret = -EINVAL; 
245 	size_t items; 
246 	int c; 
247 
 
248 	fprintf(stdout, 
249 			"\nBAT analysed signal has %d frames at %d Hz, %d channels, " 
250 			"%d bytes per sample\n", bat->frames, bat->rate, 
251 			bat->channels, bat->sample_size); 
252 
 
253 	bat->fp = fopen(bat->capture.file, "rb"); 
254 	if (bat->fp == NULL) { 
255 		loge(E_OPENFILEC, "%s", bat->capture.file); 
256 		ret = -ENOENT; 
257 		goto exit; 
258 	} 
259 
 
260 	bat->buf = malloc(bat->frames * bat->frame_size); 
261 	if (bat->buf == NULL) { 
262 		ret = -ENOMEM; 
263 		goto exit; 
264 	} 
265 
 
266 	/* Skip header */ 
267 	ret = read_wav_header(bat, bat->capture.file, bat->fp, true); 
268 	if (ret != 0) 
269 		goto exit; 
270 
 
271 	items = fread(bat->buf, bat->frame_size, bat->frames, bat->fp); 
272 	if (items != bat->frames) { 
273 		ret = -EIO; 
274 		goto exit; 
275 	} 
276 
 
277 	ret = reorder_data(bat); 
278 	if (ret != 0) 
279 		goto exit; 
280 
 
281 	for (c = 0; c < bat->channels; c++) { 
282 		struct analyze a; 
283 		printf("\nChannel %i - ", c + 1); 
284 		printf("Checking for target frequency %2.2f Hz\n", 
285 				bat->target_freq[c]); 
286 		a.buf = bat->buf + 
287 				c * bat->frames * bat->frame_size 
288 				/ bat->channels; 
289 		ret = find_and_check_harmonics(bat, &a, c); 
290 	} 
291 exit: 
292 	if (bat->buf) 
293 		free(bat->buf); 
294 	if (bat->fp) 
295 		fclose(bat->fp); 
296 	return ret; 
297 } 
