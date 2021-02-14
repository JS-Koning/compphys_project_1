# Weekly progress journal
[[_TOC_]]
## Instructions

In this journal you will document your progress of the project, making use of the weekly milestones.

Every week you should 

1. write down **on the day of the lecture** a short plan (bullet list is sufficient) of how you want to 
   reach the weekly milestones. Think about how to distribute work in the group, 
   what pieces of code functionality need to be implemented.
2. write about your progress **until Sunday, 23:59** before the next lecture with respect to the milestones.
   Substantiate your progress with links to code, pictures or test results. Reflect on the
   relation to your original plan.

We will give feedback on your progress on Tuesday before the following lecture. Consult the 
[grading scheme](https://computationalphysics.quantumtinkerer.tudelft.nl/proj1-moldyn-grading/) 
for details how the journal enters your grade.

Note that the file format of the journal is *markdown*. This is a flexible and easy method of 
converting text to HTML. 
Documentation of the syntax of markdown can be found 
[here](https://docs.gitlab.com/ee/user/markdown.html#gfm-extends-standard-markdown). 
You will find how to include [links](https://docs.gitlab.com/ee/user/markdown.html#links) and 
[images](https://docs.gitlab.com/ee/user/markdown.html#images) particularly
useful.

## Week 1
1. General
    - Get used to gitlab, make first commits and verify that these work.
    - Discuss working method, match schedules to work together.
    
2. Objectives
    - [x] Calculate the force on each particle using the Lennard-Jones potential
        - still needs checking whether this is correct with literature, mistakes might be in mass used for example.
    - [x] Implement the Euler method for time evolution
        - basic concept is implemented, however lacks basic checking due to time constraints. It runs, but might be incorrect.
    - [ ] Implement the periodic boundary condition
        - hardest part of this week, and did not know where to start. Especially since the force changes direction when a particle leaves the box. duplicating boxes seems to be an option, but then we should approach an infinite amount of boxes for the simulation to be correct? That seems incorrect.
    - [x] Define a function that calculates the total energy of the system
        - Still needs mass scaling, see pseudocode comment in the code. Would have been better implemented if we started checking our code for simulation errors, however due to boundary conditions failure this was not included. 
    - [x] Make sure that you commit your code regularly to Gitlab.
3. Things that need improvement
    - Velocity initial vector needs some improvements, this looks not right.  Also, it is based on a 3D function.
    - Simulation is still not finished totally
    - Overall code "running". Currently, it is function after function to tackle the problems, not yet a working running program. Still a lot of work in progress thus. 
4. Things that went right
    - We implemented the code for many particles, not just 2. debugging will be harder in future.
5. Things that went wrong
    - Underestimated the time requirement for programming. There is so much more to do and improve, although our code has decent functionality, it needs a whole lot more attention to details.
## Week 2
(due 21 February 2021, 23:59)


## Week 3
(due 28 February 2021, 23:59)


## Week 4
(due 7 March 2021, 23:59)


## Week 5
(due 14 March 2021, 23:59)
