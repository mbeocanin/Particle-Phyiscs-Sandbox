# MLCP Particle Phyiscs Sandbox
This is MATLAB implementation of a 2D particle physics simulation based on solving MLCP (Mixed Linear Complementarity Problem) problem using Projected Gauss-Seidel method.<br />
It is based on a model designed by Erin Catto and presented in this article: https://box2d.org/files/ErinCatto_IterativeDynamics_GDC2005.pdf

# Motivation
This project was made with the idea in mind to explore the original (pure math-based) approach to solving the MLCP problem (mentioned in the article above). It has since been surpassed by Sequential Impulse Method (https://box2d.org/files/ErinCatto_SequentialImpulses_GDC2006.pdf) which is more programmer friendly, reduces some serious performance drawbacks and provides additional flexibility in controlling the simulation compared to the original approach.<br />
Examples in this project are provided for learning purposes, since to my knowledge no other simple and compact implementation of the original approach exists elsewhere.

Since the goal was to keep the code profile as simple as possible while demonstrating the MLCP solution method, the simulation is limited to particle physics, meaning the bodies symmetrical circular shapes with no orientation and rotation.

# Running
Implementation was made and tested using MATLAB 2016a x64.<br />
To run:
1. choose src folder as MATLAB workspace directory
2. choose and run one example script at a time

GUI parameters are described through inline comments at the start of each of the example scripts.
