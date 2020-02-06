# Weekly progress journal

This is your progress journal: a simple planner that will help you make meaningful progress on your most important work.

You will use this journal to keep track of your progress in the project on an (at least) weekly basis. 
It is to be used as a reference for yourself and your team, to help you track what you are currently struggling with or what you have achieved this week.
In addition, it will also be read by the course team to review and grade your progress on the weekly milestones, listed in the lecture notes 
and in the issues that are automatically opened for you every week. Keeping a good journal is thus a mandatory part of the course!

As a guideline, a basic journal should explain what the project achieved or show a problem with the project, and show data to confirm or illustrate this.
A good journal should also outline a good plan of action for the next time.

Note that the file format of the journal is *markdown*. This is a flexible and easy method of converting text to HTML, 
and therefore very suitable for usage on a website like this. This file provides you with some basic examples of how to format your text using markdown,
but a more extensive list of options can for example be found [here](https://github.com/adam-p/markdown-here/wiki/Markdown-Cheatsheet).

## Week 1
This is a brief example of what your journal entry for week 1 could be. 

### 12-02-2020
This week we started working on project 1. We managed to form a team and started writing our first piece of code. We begun by setting up 
a structure that will store the data for our particles.
We have started writing the function that calculates the force on the particles, but we haven't managed to finish this today. 
We will come back to this later this week.
We did however make our first commit using git!

### 14-02-2020
Today we managed to implement a number of things. We wrote code that
1  implements the force on each particle using the Lennard-Jones potential
2  implements Euler's method of integration
3  wrote the code for periodic boundary conditions

The following figure illustrates the implementation of these points:
** insert figure; maybe we want to provide one as an example of how to input an image in markdown **

As seen, we plot the trajectory of two particles as a function of time. They don't seem to end up with infinitely high velocities, so that is nice!
We have not yet managed to plot the total energy of the system; we'll return to that later this week.

### 17-02-2020
Today we did some work on the formatting of our code, making it a bit easier to read. We also forgot to commit last time, which we did now. 
From here on we'll remember to do so regularly! On the bright side, we managed to calculate the total energy of the system. The figure below shows how this evolves in time:
** insert figure**

We managed to complete all of the milestones for this week, so we're happy about that! 
Our plan for the next week is to work on the concepts introduced in the next lecture.