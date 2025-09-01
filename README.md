# Cusp

Hardware-accelerated fractional Fourier transform library for Swift.

## Installation

### Swift Package Manager

Simply add Cusp to your `Package.swift` file: 

```swift
let package = Package(
    name: "Example",
    dependencies: [
        .package(url: "https://github.com/dclelland/Cusp.git", from: "0.1.0"),
    ],
    targets: [
        .target(name: "Example", dependencies: ["Cusp"])
    ]
)
```

Then import Cusp into your Swift files:

```swift
import Cusp
```

## Links

### Dependencies

- [apple/swift-numerics](https://github.com/apple/swift-numerics)
- [dclelland/Plinth](https://github.com/dclelland/Plinth)

### References/prior art

- [Digital computation of the fractional Fourier transform](https://ieeexplore.ieee.org/document/536672/) (Ozaktas et. al, 1996)
- [The Discrete Fractional Fourier Transform](https://ieeexplore.ieee.org/document/756324/) (Candan et. al, 1999)
- [tunakasif/torch_frft](https://github.com/tunakasif/torch-frft)

## Todo

- [ ] Figure out how to handle passing through the correct-sized `FFTSetup`
- [ ] Write faster 1D FrFT implementation
    - Look at caching the chirp matrices
- [ ] Write faster 2D FrFT implementation
- [ ] Fix the DFrFT Hamiltonian calculation for lengths which aren't powers of two

# Documentation

Work in progress...
