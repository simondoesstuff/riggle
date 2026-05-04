## Chuckle: World-Scale Interval Search

Simon Walker, Ryan Layer

May 4, 2026

### **Background**

Biological research relies heavily on interval sets to define drug targets, disease variants, and binding sites. To discover meaningful relationships, researchers must compute overlaps across massive databases of previously published experiments. While legacy tools like Giggle and IGD paved the way, the sheer volume of modern genomics data requires a new approach.

### **Problem Statement**

Existing interval search engines bottleneck because **they are not parallelized**. They rely on single-threaded execution and memory-heavy interval trees. As datasets scale into the billions, these legacy tools suffer from catastrophic memory spikes and exponential slowdowns. A new engine is needed to handle modern clinical scale.

### **Data Profile & Objectives**

- **The Input:** Batches of simple genomic interval files: (chr, start, end).
- **The Output Objective:** Chuckle is designed strictly for **file-to-file statistics**. Rather than extracting millions of individual coordinate overlaps, Chuckle generates a dense statistical matrix that reveals the _set-level significance_ of your query against the entire database.

---

### **Approach and Methodology**

To achieve world-scale speed, Chuckle introduces three highly optimized conceptual shifts:

**1. Shifting Complexity: The Index is Just a List**

Traditional engines waste massive amounts of memory building complex, hierarchical trees. Chuckle flips this: our index is simply a flat, sorted list of coordinates. Moving the complexity from the index phase to the query phase allows us to achieve unprecedented read speeds.

**2. The Sweep Principle (No Point Queries)**

Legacy tools treat every query interval as a blind "point" search, scanning the database from scratch every time. Chuckle recognizes that both the query and the database are pre-sorted lists. By sweeping a single line across the chromosome, we instantly tally overlaps and skip vast amounts of empty biological space without ever looking backward.

**3. True Parallelization**

Because our data is stored as a flat list rather than a tangled tree, it is easily divided. Chuckle chops the query workload into blocks and distributes them across idle processor cores. Multiple sections of the genome are swept simultaneously, merging into a final statistical matrix at the end.

---

### **Conclusion**

Chuckle successfully executes file-to-file interval search significantly faster than all benchmarked predecessors, including IGD and Giggle. Because it is fully parallelized and avoids memory bloat, it scales effortlessly into the tens of billions of datapoints.

By rethinking how interval data is stored and traversed, Chuckle provides researchers with a highly scalable, statistically driven search engine capable of mapping discoveries against the entirety of the world's experimental data on demand.

### References

    GIGGLE Layer, R. M., Pedersen, B. S., DiSera, T., Marth, G. T., Gertz, J., & Quinlan, A. R. (2018). GIGGLE: a search engine for large-scale integrated genome analysis. Nature Methods, 15(2), 123-126.[ https://doi.org/10.1038/nmeth.4556](https://doi.org/10.1038/nmeth.4556)


    IGD: high-performance search for large-scale genomic interval datasets. Bioinformatics, 37(1), 118-120.[ https://doi.org/10.1093/bioinformatics/btaa1062](https://doi.org/10.1093/bioinformatics/btaa1062)


    BEDTools: a flexible suite of utilities for comparing genomic features. Bioinformatics, 26(6), 841-842.[ https://doi.org/10.1093/bioinformatics/btq033](https://doi.org/10.1093/bioinformatics/btq033)


    Nested Containment List (NCList): a new algorithm for accelerating interval query of genome alignment and interval databases. Bioinformatics, 23(11), 1386-1393.[ https://doi.org/10.1093/bioinformatics/btl647](https://doi.org/10.1093/bioinformatics/btl647)


    AIList (Augmented Interval List) Feng, J., Ratan, A., & Sheffield, N. C. (2019). Augmented Interval List: a novel data structure for efficient genomic interval search. Bioinformatics, 35(23), 4907-4911.[ https://doi.org/10.1093/bioinformatics/btz407](https://doi.org/10.1093/bioinformatics/btz407)
