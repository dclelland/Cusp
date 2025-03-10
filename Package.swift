// swift-tools-version: 5.8
// The swift-tools-version declares the minimum version of Swift required to build this package.

import PackageDescription

let package = Package(
    name: "Cusp",
    platforms: [
        .macOS(.v13),
        .iOS(.v14)
    ],
    products: [
        .library(name: "Cusp", targets: ["Cusp"])
    ],
    dependencies: [
        .package(url: "https://github.com/dclelland/Plinth", from: "2.10.2")
    ],
    targets: [
        .target(
            name: "Cusp",
            dependencies: [
                .product(name: "Plinth", package: "Plinth")
            ]
        )
    ]
)
