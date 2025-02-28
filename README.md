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

## Todo

- [x] Add 1D FFT to Plinth
- [x] Write 1D FrFT implementation
- [x] Write 2D FrFT implementation
- [x] Write 1D LCT implementation
- [x] Write 2D LCT implementation
- [x] Write LCT matrix generators
    - Does this need arithmetic operations? Perhaps there is a generic 2x2 complex matrix type somewhere...
    
- [ ] Figure out why the FrFT and LCT are not doing the same thing; check against canonical implementations
- [ ] What is the correct way to reconstruct the original image?
- [ ] Add check in LCT for `b == 0`, as this causes a division by zero
- [ ] Should the LCT support complex-valued matrices?
- [ ] Sanity check the fractional order of 0.5 thing

- [ ] Add `Float` and `Double` implementations

# Documentation

Work in progress...
