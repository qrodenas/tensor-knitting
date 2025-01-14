# Next Goals
### 10/12 Final meeting

- Main problem now is recovering the results in qiskit, because the function reconstruct_expectation_values returns the exp += coeff*subcircuit_expval.
Options are: changing the qiskit function (should be easy), or creating my own function (perharps a bit harder but should just be a few sums)

- After that I have to think what is the best approach to integrate the qiskit and julia code for doing time evolution (maybe calling the julia code from python?)

- Look at picture in phone to know what you have to implement (comparing growth of link dimension between the approaches, and checking that both approaches yield the same expvals for time evolution).
Also some check of the runtime could be cool but I have doubts on how to implement it

- Start writing!! Inspiration from Gian's papers and the tensor network papers.  

### 26/11

- Add optimization when transpiling

- Check behaviour of MPS simulators when changing coupling constant with distance. Also check what happens if we
change link dimension.

- Fix the code for computing the expectation value with MPS (more general observable definition)

- See if there is a problem when we sample subcircuits (QPD not summing up to one)? Qiskit problem?s

- Add maxdim and cutoff to all apply().

- I am finding that maybe manually contracting tensors (manually computing inner products) can be faster

- Ask Gian where my code can be made more efficient.

- Look at older things
### 15/11

- Let’s see how big we can go on the MPO.  
  *(Remember the system looks like a kite, we want to simulate only the tail of the kite.)*
  Use MPS simulator in qiskit

- See how bond dimension in the tail works depending on JXX (constant J).

- Look at varying the coupling constant with distance *(picture)*.


## Old (mostly done)

- Try to do mixed MPS-MPO simulation, filter when there is no measurement, and only use MPO.

- **Add channels in the `build_mpo_sequence` function.**  
  Make sure you sum them!

- **Index thing in Qiskit.**  
  (Getting indices on which gates act in circuits)

- **Come up with some tests.**

- **Get expectation values for all subcircuits and compare with the full circuit run in Qiskit.**

- **Fix bug:** Y, X gates in Hamiltonian.

- **Fix poor observable definition.**  
  Good for now.

- **Make code cuter.**

- **Generalize so that it works for knitting of multiple gates.**

---

Now that we have the framework, it would be interesting to see what circuits we can simulate. For this, I would again take a Trotter circuit and try to understand:

- How big can you make the MPO circuit to still run it efficiently with good accuracy?  
  *(You can compare with the Qiskit MPS simulator if you want.)*

- Try some small bond dimensions to see how the accuracy changes with the bond dimension.

- If you fix the Hamiltonian, how many Trotter steps can you run?  
  Both just with MPO, but also with cutting a gate in every step.

---

The code works but it doesn’t seem to benefit much from using tensor networks. It's faster for more qubits and only 
one cut but for multiple cuts seems slower, at least for a few qubits. 
I think the problem is using MPOs from the start.  
We should try to avoid MPOs as much as possible and then use them at the end to introduce the measurements.

I think it’s probably time for some literature review.  
How to run quantum circuits efficiently, trying to find if someone has looked at mid-circuit measurements.  
Another way to do it is to hard-code the measurement and work always with MPS.


