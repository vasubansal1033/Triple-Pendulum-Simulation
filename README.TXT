Triple Pendulum.iynb is the Jupyter Notebook. To open that file you have to first install Jupyter Notebook
Triple Pendulum code.pdf is the pdf version of the Jupyter Notebook.

The code asks the user lengths of the pendulums and the mass at the bobs for each pendulum and the duration for which the simulation is required. Then the initial angles with the vertical and the initial angular velocity of the pendulums is asked. Based on these values and initial conditions, the code exports a video of the required duration in .mp4 format.

The code can also plot the locus of the three bobs, as seen in the image and the video.

The simulation of a Double Pendulum can also be obtained by putting M3 and L3 to be zero.

TriplePendulum_Vertical.mp4 is the simulation when the pendulum is suspended vertically upwards. 
Triple_Pendulum_Locus.mp4 demonstrates the release of pendulum horizontally, with the scatter plot of each bob.

We have also demonstrated the Butterfly Effect (read more at https://en.wikipedia.org/wiki/Butterfly_effect). The video was created by composition of three separate videos. Then a free software 'OpenShot' was used to composite them. The pendulums had initial angles of (90,90,90), (89.5,89.5,89.5) and (90.5,90.5,90.5). As we can see they show similar motion intially but it changes drastically after sometime. This due to the large dependence on the inital conditions in a chaotic system. Read more at https://simple.wikipedia.org/wiki/Chaos_theory

code.py file is the Python file which contains the code. Please make sure you have installed the required and latest versions of the libraries on your system. The code is written in Python 3. So it may not work for Python 2.x versions.
