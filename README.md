# R-BMM: A Reconfigurable Barrett Modular Multiplier Architecture for High-Performance Cryptographic Systems

Modular multiplication is a key building block in cryptographic systems, digital signatures, and secure communication, where speed, power, and hardware efficiency are all equally important. Traditional multiplier architectures are usually tailored to a fixed operand size, which not only restricts flexibility but also wastes resources when multiple bit-widths are
supported. To overcome this limitation, we propose a Recon-
figurable Barrett modular multiplier (R-BMM) that adaptively
performs one 256-bit, three 128-bit, or nine 64-bit modular
multiplications in parallel. Instead of relying on several dedicated
multipliers for different configurations, the design makes use
of a fixed set of reusable multiplier core, significantly cutting
down on hardware overhead while maintaining superior perfor-
mance. The proposed architecture is implemented, synthesized
and evaluated on both FPGA ASIC to validate its efficiency.
On FPGA, the design achieves a 52.178% reduction in area,
an impressive 57.33% reduction in power consumption, and a
49.142% reduction in critical path delay compared to state-of-the-
art designs. ASIC synthesis results obtained using Cadence Genus
tool confirm these trends, demonstrating substantial savings in
gate count and dynamic power, together with a higher achievable
frequency, thereby reinforcing the scalability of the proposed
approach. These consistent improvements across FPGA and
ASIC platforms makes the architecture particularly attractive for
resource-constrained cryptographic accelerators and embedded
platforms, where efficiency and adaptability are crucial. Overall,
the proposed reconfigurable multiplier offers a hardware efficient
and high-performance solution for future secure computing
applications. The designs are made freely available for easy
adoption and further usage to the researchers and designers
community.
