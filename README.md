# SystemC Test Bench for Canny Edge Detector

This project demonstrates a SystemC implementation of a Canny Edge Detector with a hierarchical test bench architecture for embedded systems modeling and design.

## Environment

Linux environment on UC Irvine's remote server.

## Architecture

The system employs a structured test bench approach with clear module separation:

### Top-Level Module Hierarchy
```
Top
├── Stimulus (input file reader)
├── Platform (processing platform)
│   ├── DataIn (input interface)
│   ├── DUT (Canny algorithm core)
│   └── DataOut (output interface)
└── Monitor (output file writer)
```

### Communication Infrastructure
- **FIFO Channels**: Uses `sc_fifo<IMAGE>` for inter-module communication
- **Custom Image Type**: Implements IMAGE struct with proper copy operators for array handling
- **Buffer Size**: Single-slot FIFO queues (size=1) for minimal latency

## SystemC Features Demonstrated

### Module Design Patterns
- **SC_MODULE** declarations for all components
- **SC_THREAD** processes for concurrent execution
- Proper **SC_CTOR** constructors with stack size management

### Communication Mechanisms
- **sc_fifo_in/sc_fifo_out** ports for type-safe communication
- Port binding in `before_end_of_elaboration()` phase
- Hierarchical channel instantiation

## Technical Implementation

### Image Processing
- **Resolution**: 2704×1520 pixels
- **Format**: PGM (Portable Graymap)
- **Processing**: 30-frame video sequence
- **Memory**: Static array allocation with custom wrapper

### Stack Size Configuration

Large images stored in local variables require increasing stack size for SystemC models:

1. **Root thread stack size**
   - For `csh` / `tcsh`:
     ```sh
     limit stacksize 128 megabytes
     ```
   - For `sh` / `bash`:
     ```sh
     ulimit -s 128000
     ```

2. **SC_THREAD stack size**  
   After each `SC_THREAD()` declaration, add:
   ```cpp
   set_stack_size(128*1024*1024);

## Running the Project

### Build and Execute
```bash
# Compile the project
make

# Run the Canny edge detector
./canny

# Clean build artifacts
make clean
```

### Test with Sample Data
```bash
# Ensure input images are in the correct directory
# Run the complete processing pipeline
make run

# Compare output with reference images (if available)
make test
```