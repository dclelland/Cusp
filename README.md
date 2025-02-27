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

- [ ] Add 1D FFT to Plinth
- [ ] Write 1D FrFT implementation
- [ ] Write 2D FrFT implementation
- [ ] Write 1D LCT implementation
- [ ] Write 2D LCT implementation
- [ ] Write LCT matrix generators
    - Does this need arithmetic operations? Perhaps there is a generic 2x2 complex matrix type somewhere...

# Documentation

Work in progress...
