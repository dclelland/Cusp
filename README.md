# Cusp

Hardware-accelerated linear canonical transformation library for Swift.

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

- [ ] What is the correct way to reconstruct the original image?
- [ ] Add check in LCT for `b == 0`, as this causes a division by zero
- [ ] Should the LCT support complex-valued matrices?
- [ ] Sanity check the fractional order of 0.5 thing

- [ ] Add `Float` and `Double` implementations

# Documentation

Work in progress...
