# Parallelization_Energy_Calculation_water_molecules
<<<<<<< HEAD
This project was based on parallelizing the calculation of 10k&amp;100k water molecules with different box sizes. The results showed some consistency except for few irregularities in overhead and other fields, but all in all, the results were a good practical example to reflect the time efficiency of parallelization.
=======
This project was based on parallelizing the calculation of 10k & 100k water molecules with different box sizes. The results showed some consistency except for few irregularities in overhead and other fields, but all in all, the results were a good practical example to reflect the time efficiency of parallelization.

# Methods:

## 1. MPI method:



### Reading file

As for reading of files, Sequential reading is reliable for any number of threads with small changes compared to serial code. However, in this mini project, we are comparing the timings of the three different methods which shows how option 3 and option 1 are quite similar since they aren’t affected by the number of cores. On the other hand, sequential reading time of option 2 is proportional to the number of cores, therefore it’s very slow when going for higher number of cores.

![](reading.png)



### Timing and efficiency


For parallel calculation, numbers look relatively and steadily improving compared to the number of cores. There are – of course – some unexpected irregularities such as the code of 10k.gro with 24 cores on a single machine.


![](timing.png)

### Graph analysis

It had an overhead of 8 compared to 1 in 4 cores, I further inspected the issue and run few “qsub” only for that specific number of cores, but I still get the same result. With the single core code in parallel and serial, there is a tiny difference which is expected due to load and some extra printing tasks in the parallel one.

![](graph1.png)



For the imbalance, it wasn’t zero all the time except for the 1 core, but it’s extremely slight error, and it increases linearly with the more cores added (as shown in the graph).


![](graph2.png)

most of the parallelized codes worked well, the outputs were mostly linearly proportional and close to the theoretical time itself, proving that parallelization is a very important skill to efficiently accomplish heavy workload.
>>>>>>> a1bcab3 (adding MPI files)
