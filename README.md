# SOFA-object-for-Pure-data
SOFA-object-for-Pure-data is a graphical programming language called Pure-data that uses a Spatial Oriented Format (SOFA) file to convolve a user-specified .wav file with IR (impulse response) based on elevation azimuth and distance position information.
It is an object that performs IR (impulse response) convolution according to the user-specified .wav file and elevation azimuth and distance position information, and outputs the result.
With this object, you can manipulate stereophonic sound in real time within Pure-data.


# Requirement

* Pure-data
* fftw-3.3.9
* libmysofa



# Installation
Pure-data
http://puredata.info/downloads/pure-data
The site to install Puredata.
Please install Puredata according to your OS here.

fftw
http://www.fftw.org/download.html
Install the latest fftw here

> $ tar -zxvf fftw-3.3.9.tar.gz

> $ cd fftw-3.3.9

> $ make

> $ sudo make install

Installation is complete.

Make sure that libfftw3.a, libfftw3.la, libfftw3f.a, libfftw3f.la, and libmysofa.a are in ~/usr/local/lib/.

If you don't have them, move them to ~/usr/local/lib.


libmysofa
https://github.com/hoene/libmysofa
Follow the Compile here.

If you do not have Cunit installed, please install it.
https://sourceforge.net/projects/cunit/files/

Also, if you do not have Node.js installed, please install it.


# Usage
Execute this command
Put the made file in the path of pure-data.
Compilation is complete.

How to use pdsofa~help.pd

> $ cd mysofa~

> $ make

First, press the bung connected to the [openpanel] and select the sofafile you want to include.

Next, select the .wav file you want to put in the sound file array.

Turn on the DSP.

Press the bang connected to [tabplay~].

Change the azimuth, elevation, and distance values, and the sound will be output accordingly.

There is a video for reference.

https://youtu.be/DrBP1Iz9Ots 
(English)

https://youtu.be/eYMlJGgj1rw 
(Japan)

# Author
* Author Sakuya Fujisawa
* E-mail m5251103@u-aizu.ac.jp

