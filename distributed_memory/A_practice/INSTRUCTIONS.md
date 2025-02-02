# Exercise: Estimating Pi Using Numerical Integration with MPI

## Objective

In this exercise, you will estimate the value of $\pi$ using numerical integration and parallel processing with MPI. The method involves approximating the integral of the function:

$$
f(x) = \frac{4}{1 + x^2}
$$

over the interval $[0, 1]$, which is equal to $\pi$.

---

## Instructions

1. **Divide the Domain:**
   - Split the interval $[0, 1]$ into $N$ subintervals (e.g., $N = 10^6$).
   - Each MPI process should handle a portion of the subintervals.

2. **Compute Local Contributions:**
   - Each process computes the partial sum of the function values for its assigned subintervals.

3. **Communicate and Aggregate Results:**
   - Use MPI functions to collect partial sums from all processes and compute the final approximation of π on the root process.

4. **Implement the Program:**
   - Write the MPI program using `MPI_Init`, `MPI_Comm_size`, `MPI_Comm_rank`, and `MPI_Reduce`.
   - Use a parallel numerical integration technique such as the midpoint rule or trapezoidal rule.

5. **Test and Experiment:**
   - Test your implementation with different values of $N$ and measure the time taken for computation.
   - Observe how the number of MPI processes affects performance.

6. **Output:**
   - Display the estimated value of $\pi$, the actual error ($|\pi - \text{estimated value}|$), and the computation time for each test.

---

## Example Output

For $N = 10^6$ and 4 MPI processes:

```Estimated value of pi: 3.14159265 Actual error: 0.00000001 Time taken: 0.05 seconds```


# Exercise: Parallel Matrix - Vector Multiplication Over Rows Using MPI

## Objective

In this exercise, you will implement matrix -vector multiplication in parallel using MPI. The computation will be distributed across the rows of the first matrix ($A$), with each MPI process computing a subset of the rows in the resulting vector ($\mathbf{A} x = b$).

---

## Problem Description

Given:
- $\mathbf{A}$: An $N \times K$ matrix
- $x$: A $K \times 1$ vector

The result is vector $b$, an $N \times 1$ matrix, where:

$$
b_i = \sum_{k=0}^{K-1} A_{ik}x_k
$$

---

## Instructions

1. **Initialize MPI:**
   - Use `MPI_Init`, `MPI_Comm_size`, and `MPI_Comm_rank` to set up the MPI environment.
   - Determine the total number of MPI processes ($P$).

2. **Divide Rows Among Processes:**
   - Split the rows of matrix $\mathbf{A}$ evenly among the $P$ processes. 
   - Each process handles $N / P$ rows (or more, for uneven divisions).

3. **Broadcast Vector $x$:**
   - Vector $x$ is required by all processes for the computation.
   - Use `MPI_Bcast` or initialize the same vector on every process to ensure $x$ is in all processes.

4. **Compute Local Contributions:**
   - Each process calculates the rows of $b$ assigned to it using its subset of rows from $\mathbf{A}$ and the entire $x$ vector.

5. **Gather Results:**
   - Use `MPI_Gather` (or `MPI_Gatherv` for uneven distributions) to collect the computed rows of $b$ from all processes on the root process.

6. **Finalize MPI:**
   - Use `MPI_Finalize` to clean up.

---

## Implementation Details

- **Input:**
  - Generate or read matrix $\mathbf{A}$ and vector $x$.
  - Ensure the dimensions of $\mathbf{A}$ and $x$ are compatible for multiplication.

- **Output:**
  - Display or save the resulting vector $b$ on the root process.
  - Print the computation time for performance analysis.

- **Optimization:**
  - Consider minimizing communication overhead by dividing rows more efficiently.
  - Use collective communication functions where possible.
---

# Exercise: Ring Traversal Using MPI

## Objective

In this exercise, you will implement a ring traversal using MPI. Each MPI process will send its rank value to the next process in the ring, and the message will circulate until it returns to the original process.

---

## Problem Description

Given $P$ MPI processes, arrange them in a logical ring. Each process will:
1. Send its rank value to the next process in the ring.
2. Receive a rank value from the previous process.
3. Continue forwarding the rank value until the message completes a full circle and returns to the originating process.

---

## Instructions

1. **Initialize MPI:**
   - Use `MPI_Init`, `MPI_Comm_size`, and `MPI_Comm_rank` to set up the MPI environment.
   - Determine the total number of MPI processes ($P$).

2. **Define the Ring:**
   - Each process should send data to `(rank + 1) % P`.
   - Each process should receive data from `(rank - 1 + P) % P` to handle the circular structure.

3. **Implement the Message Passing:**
   - Use `MPI_Send` and `MPI_Recv` or `MPI_Sendrecv` to send and receive rank values.
   - Ensure proper synchronization by setting up matching send and receive operations.

4. **Handle the Full Circle:**
   - Each process forwards the received rank value to the next process.
   - When the message returns to the originating process, print the full traversal order.

5. **Finalize MPI:**
   - Use `MPI_Finalize` to clean up.

---

## Example Execution

For $P = 4$ MPI processes:

```
Process 0 received rank 3.
Process 1 received rank 0.
Process 2 received rank 1.
Process 3 received rank 2.
```

## Extra Step
Try to implement this using a ```MPI_Cart_comm``` communicator.

# Exercise: Estimating π Using Monte Carlo Simulation in MPI

## Objective

Estimate the value of π using a Monte Carlo simulation. Each MPI process generates random points, counts how many fall within a unit circle, and combines the results to compute π.

---

## Problem Description

We approximate the value of π by generating random points within a square and counting how many fall inside the unit circle. The relationship is:
$$
\pi \approx 4 \times \frac{\text{Points inside the circle}}{\text{Total points}}
$$
Using $P$ MPI processes:
1. Divide the total number of points $N$ evenly across processes.
2. Each process generates random points and counts how many are inside the circle.
3. Use a collective operation to sum the counts and calculate π.

---

## Instructions

1. **Initialize MPI:**
   - Set up the MPI environment using `MPI_Init`, `MPI_Comm_size`, and `MPI_Comm_rank`.
   - Determine the number of processes ($P$) and the rank of each process.

2. **Distribute the Work:**
   - Choose the total number of points ($N$).
   - Divide $N$ evenly among all processes, so each process computes $N/P$ points.

3. **Generate Random Points:**
   - Each process generates random $(x, y)$ coordinates within the range $[-1, 1]$.
   - Use the condition $x^2 + y^2 \leq 1$ to determine if a point lies inside the unit circle.

4. **Count Points:**
   - Each process keeps a local count of points inside the circle.

5. **Reduce Results:**
   - Use `MPI_Reduce` to sum up the counts from all processes at the root process.
   - The root process computes the final value of π.

6. **Print Results:**
   - Each process can print its local results.
   - The root process prints the estimated value of π and the actual error compared to the true value.

7. **Finalize MPI:**
   - Use `MPI_Finalize` to clean up the MPI environment.

---

## Example Execution

For $N = 10^6$ and $P = 4$ processes:

- Each process computes $N/P = 250,000$ points.
- Output might look like this:

```
Process 0: Local points = 250000, Inside circle = 196485 
Process 1: Local points = 250000, Inside circle = 196321 
Process 2: Local points = 250000, Inside circle = 196401 
Process 3: Local points = 250000, Inside circle = 196392
Total points = 1000000 Estimated value of π = 3.141504 Actual error = 0.000089
```

---

## Implementation Details

- **Input:**
  - Total number of points ($N$).
  - Number of processes ($P$).

- **Output:**
  - Each process prints its local results (points generated, points inside circle).
  - The root process prints the final estimated value of π.

- **Edge Cases:**
  - $P = 1$: All points are computed by a single process.
  - Large $N$: Ensure random number generation remains efficient.

---

## Extensions

1. **Measure Performance:**
   - Use `MPI_Wtime` to measure the total computation time.
   - Compare performance as $P$ increases.

2. **Improve Accuracy:**
   - Experiment with larger values of $N$.
   - Use a higher-quality random number generator.

3. **Analyze Scaling:**
   - Run the simulation with varying numbers of processes and observe scaling behavior.

4. **Visualize Results:**
   - Output a small subset of points to visualize their distribution.

---