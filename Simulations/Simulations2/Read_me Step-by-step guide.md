# Replicating Bivariate Benchmark Comparison Results

This guide explains how to reproduce the results from the **Bivariate Benchmark Comparison** section of the manuscript.

> ⚠️ I’m more of a mathematician than a computer scientist, so Python–R fusion here is more “duct tape” than “art”; but it runs on my machine. Email me if yours doesnt work: juraj.bodik@unil.ch
---

## Steps

1. **Upload the CPCM function**

   * Open `CPCM_function.R` and run the file (this loads the main CPCM function).

2. **Download required repositories**

   * [bQCD-master.zip](https://github.com/tagas/bQCD)
   * [loci-master.zip](https://github.com/AlexImmer/loci)

3. **Run Python code for the methods**

   * Extract `loci-master.zip`.
   * Open a Python session/module.
   * Set the working directory to the extracted `loci-master` folder.
   * Run the following code:

     ```python
     from causa.datasets import AN, LS, MNU, SIMG, ANs, CausalDataset, Tuebingen, SIM, LSs
     from causa.heci import HECI
     from causa.loci import loci
     ```

4. **Run R code for the methods**

   * Extract `bQCD-master.zip`.
   * Set your R working directory to the extracted `bQCD-master` folder.
   * Run `Data_generators.R` and `Baseline_methods_in_R.R`.

5. **Run the simulation**

   * Open `Simulations 2.txt` and run it.

---

## Notes

* For cleaner implementations of all methods (except CPCM), use:

  * Python methods: [loci GitHub repo](https://github.com/AlexImmer/loci)
  * R methods: [bQCD GitHub repo](https://github.com/tagas/bQCD)



