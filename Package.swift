// swift-tools-version: 6.4

import PackageDescription

let package = Package(
    name: "swift-acoustics-cavity",
    platforms: [
        .macOS(.v27),
        .iOS(.v27),
        .tvOS(.v27),
        .watchOS(.v27),
        .visionOS(.v27),
    ],
    products: [
        .library(name: "Acoustics Cavity", targets: ["Acoustics Cavity"])
    ],
    dependencies: [
        .package(url: "https://github.com/swift-molecules/swift-acoustics.git", branch: "main"),
        .package(url: "https://github.com/swift-molecules/swift-galerkin.git", branch: "main"),
        .package(url: "https://github.com/swift-molecules/swift-sparse.git", branch: "main"),
        .package(url: "https://github.com/swift-molecules/swift-eigen.git", branch: "main"),
        .package(url: "https://github.com/swift-atoms/swift-point.git", branch: "main"),
        .package(url: "https://github.com/swift-atoms/swift-space.git", branch: "main"),
        .package(url: "https://github.com/swift-atoms/swift-frequency.git", branch: "main"),
        .package(url: "https://github.com/swift-atoms/swift-mass.git", branch: "main"),
        .package(url: "https://github.com/swift-molecules/swift-acoustics-plate.git", branch: "main"),
        .package(url: "https://github.com/swift-molecules/swift-acoustics-shell.git", branch: "main"),
        .package(url: "https://github.com/swift-atoms/swift-mesh.git", branch: "main"),
    ],
    targets: [
        .target(
            name: "Acoustics Cavity",
            dependencies: [
                .product(name: "Acoustics", package: "swift-acoustics"),
                .product(name: "Galerkin", package: "swift-galerkin"),
                .product(name: "Sparse", package: "swift-sparse"),
                .product(name: "Eigen", package: "swift-eigen"),
                .product(name: "Point", package: "swift-point"),
                .product(name: "Space", package: "swift-space"),
                .product(name: "Frequency", package: "swift-frequency"),
                .product(name: "Mass", package: "swift-mass"),
            ]
        ),
        .testTarget(
            name: "Acoustics Cavity Tests",
            dependencies: [
                .target(name: "Acoustics Cavity"),
                .product(name: "Acoustics Plate", package: "swift-acoustics-plate"),
                .product(name: "Acoustics Shell", package: "swift-acoustics-shell"),
                .product(name: "Mesh", package: "swift-mesh"),
            ]
        ),
    ],
    swiftLanguageModes: [.v6]
)

for target in package.targets where ![.system, .binary, .plugin, .macro].contains(target.type) {
    target.swiftSettings =
        (target.swiftSettings ?? []) + [
            .strictMemorySafety(),
            .enableUpcomingFeature("ExistentialAny"),
            .enableUpcomingFeature("InternalImportsByDefault"),
            .enableUpcomingFeature("MemberImportVisibility"),
            .enableUpcomingFeature("NonisolatedNonsendingByDefault"),
            .enableUpcomingFeature("InferIsolatedConformances"),
            .enableExperimentalFeature("Lifetimes"),
            .treatAllWarnings(as: .error),
        ]
}
