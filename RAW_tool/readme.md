#Informations

This folder contain only a first RAW version of the laboratory tool to test the back end algorithms. Their is **no graphical interface** coming with this software actually.

This code is made actually in C and may need small adjustement to be compiled in C++. Developement and tests were only made on **Linux** actually (the author work only on Linux), to compile this code on Windows some modification must be done.

#Compilation

Use the following command line to compile the code: 

    gcc -g test.c -o test -lm -lpthread

#Use

First launch the probe emulator on your PC (use the plate image for example, so the image is allways the same) and launch ./test. 

The code wait for the first line of an image and acquire only. Then it determine the envelope of each line and do the scan conversion of the image. The image is save on a txt file. 

The image can be see we octave for example with the plot_result.m file.