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

- [ ] Implement fractional Fourier transform based on [`torch_frft`](https://github.com/tunakasif/torch-frft)
    - [x] Remove horizontal shift
    - [x] Try removing the thing which multiplies the first value by 2.0 for no reason (why does it do that?)
    - [ ] Fix issue where the `reversed()` call is offset
    - [ ] Tidy up implementation
        - [ ] Extract chirp functions and compare with the LCT ones
    - [ ] 2D FFT version
    - [ ] Can you port the current changes back over to the LCT version?

# Documentation

Work in progress...
