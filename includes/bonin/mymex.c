#include <mex.h>
#include <matrix.h>
#include <stdio.h>
// [frame]=formatframe(raw,linestarts,xii,width,linelen,delay)
//delay=maxdelay-delay

void mexFunction(int nlhs, mxArray *plhs[],
		 int nrhs, const mxArray *prhs[])
{

    /* description of the variables below upon initialization */
    
    mwSize frame_dims[]={2,2};
    unsigned short i,j;
    signed short *pout, *pinCol, *poutCol;
    //unsigned char *pin;
    unsigned long *pin;
    unsigned long *pxii,*pwidth;
    unsigned long *pos, *plen, *del;
    //unsigned long *pos, *plen;
    //signed long *del;
    unsigned long count = 0;
    size_t line_n, ind_n, hf_line_n;
    signed long output_offset,input_offset;   
    
    
    /* Check for proper number of input and output arguments */    
    
    if (nrhs != 6) 
    {
	mexErrMsgTxt("Six input arguments required.");
    }
    
    if (nlhs > 1)
    {
	mexErrMsgTxt("Too many output arguments.");
    }
    
    if (mxGetClassID(prhs[0])!= mxUINT8_CLASS) 
    {
        mexErrMsgTxt("first argument should be uint8");
    }

    if (mxGetClassID(prhs[1])!= mxUINT32_CLASS)  
    {
        mexErrMsgTxt("second argument should be uint32");
    }

    if (mxGetClassID(prhs[2])!= mxUINT16_CLASS) 
    {
        mexErrMsgTxt("third argument should be uint16");
    }

    if (mxGetClassID(prhs[3])!= mxUINT16_CLASS) 
    {
        mexErrMsgTxt("Fourth argument should be uint16");
    }
    
    if (mxGetClassID(prhs[4])!= mxUINT32_CLASS) 
    {
        mexErrMsgTxt("Fifth argument should be uint32");
    }
    
    hf_line_n=mxGetNumberOfElements(prhs[1]);  // number of input cycles
    line_n = (unsigned long)mxGetNumberOfElements(prhs[1])*2;  // number of output lines
    ind_n = mxGetNumberOfElements(prhs[2]); // number of samples per line
    
   
    
    /* get inputs ptrs */
    
    
    //pin = (signed short *) mxGetPr(prhs[0]); // raw data array
    pin = (unsigned short *) mxGetPr(prhs[0]); // raw data array
    pos = (unsigned long *) mxGetPr(prhs[1]); // line start indices
    //pxii = (unsigned short *) mxGetPr(prhs[2]); // bin indices
    pxii = (unsigned long *) mxGetPr(prhs[2]); // bin indices
    pwidth = (unsigned long *) mxGetPr(prhs[3]); // output frame width
    plen = (unsigned long *) mxGetPr(prhs[4]); // input line length
    //del = (signed long *) mxGetPr(prhs[5]); // interline delay
    del = (unsigned long *) mxGetPr(prhs[5]); // interline delay

           
    /* output array line_n X *pwidth */
    
    frame_dims[0] = (mwSize) line_n; 
    frame_dims[1] = *pwidth; 
    plhs[0] = mxCreateNumericArray(2, frame_dims, mxINT16_CLASS, mxREAL);    
    //pout = (signed short *) mxGetPr(plhs[0]);
    pout = (signed short *) mxGetPr(plhs[0]);

    count = 0;
   
    
    // for 0<=j<ind_n-1
    for(j=0;j<ind_n-1;j++) // loop over columns
    { 
        
        // 1. sum samples into bins
        for(i=0;i<hf_line_n;i++) //loop over line pairs (cycles)
        {
            
            // bin samples odd lines (left-to-right part of the cycle)
            input_offset = pos[i]+*(del)+j;   
            //input_offset = (short)input_offset;
            
            if (input_offset < 0) {
                input_offset = 0;
            }
                        
            output_offset = (pxii[j]-1)*line_n+i*2; 
            
            //*(pout+output_offset) = *(pout+output_offset) + *(pin+input_offset);
            
            /*
            // bin samples even lines (right-to-left part of the cycle)
            input_offset = pos[i]+*(plen)-1+j;
            input_offset = (short)input_offset;
            
            if (input_offset < 0) {
                input_offset = 0;
            }
            output_offset = (*(pwidth)-pxii[j])*line_n+i*2+1;
            *(pout+output_offset) = *(pout+output_offset) + *(pin+input_offset);
            */
        }
        count = count + 1;

        /*
        // 2. compute average for each bin
        if (j + 1 < ind_n & pxii[j+1]>pxii[j] & count > 0) 
        {
            if (count>1)
            {
                 for(i=0;i<hf_line_n;i++)
                 {
                      *(pout+(pxii[j]-1)*line_n+i*2) = *(pout+(pxii[j]-1)*line_n+i*2) / count;
            		*(pout+(*(pwidth)-(pxii[j]-1)-1)*line_n+i*2+1) = *(pout+(*(pwidth)-(pxii[j]-1)-1)*line_n+i*2+1) / count;
                 }
            }            
             count = 0;            
        }  
         */
    }

    /*
//     // for j = ind_n-1
    j = ind_n-1;
    if (count>1)
    {
        // loop over line pairs
        for(i=0;i<hf_line_n;i++)
        {
           *(pout+(pxii[j]-1)*line_n+i*2) = *(pout+(pxii[j]-1)*line_n+i*2) / count;
           *(pout+(*(pwidth)-(pxii[j]-1)-1)*line_n+i*2+1) = *(pout+(*(pwidth)-(pxii[j]-1)-1)*line_n+i*2+1) / count;
         }
     }  
     
     */
    
    return;
}





//             fprintf("Output offset is %d",output_offset);
