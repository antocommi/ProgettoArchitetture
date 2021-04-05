# Assembly optimizations
Authors: Alessandro Barbieri, Antonella Sola, Antonio Commisso

We built a system that is able to calculate the K Nearest Neighbours of a given query x. We implemented two types of inference algorithms based on the techniques described 
in the following paper: https://hal.inria.fr/inria-00514462/document

We have two types of inference:
- Exact 
- Approximate 

The program is written in C but we optimized some part in Assembly. The algorithms use several times matrix operation so we wrote general method to deal with matrices easily 
and we drastically reduced the coding time. Furthermore we optimize the code adapting well known techniques such as: 
- Loop unrolling,
- Cache blocking,
- Etc...
