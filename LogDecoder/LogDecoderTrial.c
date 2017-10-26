/* Nicole Ortega
 * nortega1@jhu.edu
 * 8/12/2015
 *
 * Program reads a txt file and decodes each character to extract data. 
 * 8/12: Currently data is encoded in the following order: 
 * 1st byte: tap and buzz as 1st and 2nd MSB
 * 2nd - 5th byte: float time
 * 6th - 7th byte: int16 encoder count
 *
 * Ouput will be of the form,
 * buzz, tap, time, encoder counts
 * 
 * Output Example: 
 *  0, 0, 0.001, 23
 *  1, 0, 0.002, 25
 *  0, 1,0.003, 30
 * ......
 *
 * EDIT 6/30/2017: capable of decoding all files in one folder 
 * 
 * Compilation instructions: 
 *      make (or gmake. uses makefile)
 *      make clean (cleans executible)
 * Run executible: 
 *      ./LogDecoder2 date subject_num num_trials
 * OR 	LogDecoder2.exe date subject_num num_trials
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h> 
#define BIT_MASK_8 255

 /* ---------- decoding ------------ */
void chartoByte(char input, char *output)
{
    for (int i = 7; i >= 0; --i)
    {
        output[7-i] = (input & (1 << i)) ? '1' : '0' ;
      
    }
    output[8] = '\0';
}
float decodeFloat(char * input)
{
    unsigned int f = 0;

    for (int i = 3; i >= 0; i--) {
        f |= ((*(unsigned char *) input) << 8*i) & (BIT_MASK_8 << 8*i);
        input++;
    }
    
    return *((float *) &f);
}

int16_t decodeIntsixteen(char * input)
{
    unsigned int f = 0;

    for (int i = 1; i >= 0; i--) {
        f |= ((*(unsigned char *) input) << 8*i) & (BIT_MASK_8 << 8*i);
        input++;
    }
    
    return *((int16_t *) &f);
}

// A simple atoi() function from GeeksforGeeks
int myAtoi(char *str)
{
    int res = 0; // Initialize result
  
    // Iterate through all characters of input string and
    // update result
    for (int i = 0; str[i] != '\0'; ++i)
        res = res*10 + str[i] - '0';
  
    // return result.
    return res;
}

//takes in 3 arguments: 1) date 2) subject number 3) single trial

int main(int argc, char *argv[])
{   
    int trial = myAtoi(argv[3]);
	printf("\n%d\n", trial);
	char input[50];
	char output[50];
	sprintf(input, "%s/%s_Subject_%s.%d_data.txt", argv[1], argv[1], argv[2], trial);
	sprintf(output, "%s/%s_Subject_%s.%d_Res.txt", argv[1], argv[1], argv[2], trial);
	printf("%s", input);
	printf("%s\n", output);

	
    FILE *logFile = fopen(input, "rb");
    FILE *decodedFile = fopen(output, "wb");

    char tapbuzz[8];
    float time = 0;
    int16_t encoder_counts = 0;
    char buff[7];
    //unsigned char c; // Temp holder for each character in the file

    //fgetc(logFile); // Ignore the first character

    while (!feof(logFile)) 
    {
        for (int i = 0; i < 7; i++) 
        {
            buff[i] = fgetc(logFile);
        }

        chartoByte(buff[0], tapbuzz);
        time = decodeFloat(&(buff[1]));
        encoder_counts = decodeIntsixteen(&(buff[5]));
        fprintf(decodedFile, "%c %c %f %d\n", tapbuzz[0], tapbuzz[1], time, encoder_counts);
    }

    fclose(logFile);
    fclose(decodedFile);
        
    return 0;
}
