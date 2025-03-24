
# Design & Development of IEEE 802.11a Architecture on Zedboard

## Overview

This repository contains the design and development work for implementing the IEEE 802.11a wireless communication standard on a Zedboard platform. The project demonstrates the integration of high-speed digital communication systems using FPGA hardware and explores hardware/software co-design approaches to achieve robust wireless connectivity.

## Features

- **IEEE 802.11a Implementation:** Realize core functionalities of the IEEE 802.11a standard including modulation, demodulation, and signal processing.
- **FPGA-Based Design:** Developed using VHDL/Verilog for hardware modules tailored for high-speed wireless communication.
- **Embedded Software Integration:** Utilizes C/C++ code running on the ARM processor embedded in the Zedboard to manage system configuration and control.
- **Testing & Validation:** Comprehensive simulation and hardware testing to ensure protocol compliance and optimal performance.

## System Architecture

The project is divided into several key components:

- **Hardware Design:** 
  - Implements the physical (PHY) layer functionalities such as OFDM modulation/demodulation and error correction.
  - Uses Xilinx design tools (Vivado or ISE) for synthesis, implementation, and simulation.

- **Software Integration:** 
  - Contains firmware written in C/C++ for ARM-based system management.
  - Handles configuration, control, and runtime monitoring of the wireless communication modules.

- **Testing & Validation:**
  - Includes simulation test benches and on-hardware evaluation to verify design integrity.
  - Uses debugging interfaces such as serial communication for real-time data monitoring.

## Getting Started

### Prerequisites

- **Hardware:**
  - Zedboard (or a compatible Xilinx Zynq-7000 based board)
- **Software Tools:**
  - Xilinx Vivado/ISE Design Suite for FPGA development
  - Xilinx SDK for ARM development
  - Simulation tools (e.g., ModelSim)
- **Skills:**
  - Basic knowledge of digital signal processing, FPGA design, and wireless communication protocols

### Installation & Setup

1. **Clone the Repository:**

   ```bash
   git clone https://github.com/MINTU296/Design-Development-of-IEEE-802.11a-Architecture-on-Zedboard.git
   ```

2. **Open the Project:**
   - Launch the project in Xilinx Vivado or ISE using the provided project file (.xpr or .xise).

3. **Synthesize & Implement:**
   - Run synthesis and implementation steps following the standard procedures in your chosen tool.
   - Generate the FPGA bitstream after successful synthesis and implementation.

4. **Deploy on Zedboard:**
   - Use the Xilinx tools to program the Zedboard with the generated bitstream.

5. **Software Setup:**
   - Open the accompanying software project in Xilinx SDK.
   - Compile and deploy the firmware to the ARM processor on the Zedboard.

### Running the Project

- Detailed usage instructions, including simulation and debugging guidelines, are provided in the project's documentation folder.
- Monitor the system outputs via the serial interface or any configured debugging channels.

## Design Details

- **IEEE 802.11a Standard:** An overview of the standard, highlighting key protocol requirements and performance targets.
- **Module Breakdown:** Detailed descriptions of individual hardware modules such as OFDM modulators/demodulators, error correction circuits, and synchronization units.
- **Software Control:** Explanation of how the embedded software interacts with hardware modules to manage real-time operations.
- **Validation Methods:** Information on simulation test benches and hardware testing strategies used to verify the design.

## Documentation

Additional design documents, simulation results, and user manuals can be found in the [Documentation](./Documentation) folder.

## Contributors

- [MINTU296](https://github.com/MINTU296)
