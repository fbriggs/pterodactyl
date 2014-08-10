/****** NOTICE: LIMITATIONS ON USE AND DISTRIBUTION ********\
 
Ê This software is provided on an as-is basis, it may be 
Ê distributed in an unlimited fashion without charge provided
Ê that it is used for scholarly purposes - it is not for 
Ê commercial use or profit. This notice must remain unaltered. 
 
Ê Software by Dr Fred DePiero - CalPoly State University
 
\******************** END OF NOTICE ************************/ 

// header of wav file 
typedef struct{
   char rID[4];            // 'RIFF' 
   int32_t rLen;
      
   char wID[4];            // 'WAVE' 
      
   char fId[4];            // 'fmt ' 
   int32_t pcm_header_len;   // varies... 
   int16_t wFormatTag;
   int16_t nChannels;      // 1,2 for stereo data is (l,r) pairs 
   int32_t nSamplesPerSec;
   int32_t nAvgBytesPerSec;
   int16_t nBlockAlign;      
   int16_t nBitsPerSample;
}   WAV_HDR;
 
   
// header of wav file 
typedef struct{
   char dId[4];            // 'data' or 'fact' 
   int32_t dLen;
//   unsigned char *data; 
}   CHUNK_HDR;
 

// modification by Forrest Briggs 2013:
// this WAV in code will break if it is compiled as 64-bit, because it uses
// sizeof() on the above structures. the defines below are the 32-bit size of these headers

#define SIZEOF_WAV_HDR	36
#define SIZEOF_CHUNK_HDR 8

// Note: in addition to defining these sizes, I changed the definitions of WAV_HDR and CHUNK_HDR to use explicitly sized int types (i.e. int32_t and int16_t), rather 
// than long int and short int.