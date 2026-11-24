extension Acoustics.Cavity {
    /// One acoustic mode of a body of air: its frequency, its pressure shape, and the volume that sets how strongly
    /// the walls drive it.
    public struct Mode {
        /// The natural frequency.
        public let frequency: Frequency.Value<Double>
        /// The pressure at each vertex of the cavity, in its order, scaled so that the largest is one and positive.
        /// An open vertex has zero.
        public let shape: [Double]
        /// The modal volume at that scale, Λ = ∫ shape² dV. A wall moving with normal acceleration a drives the
        /// mode's amplitude q, the pressure q × shape, by q̈ + ω² q = −(ρ c² / Λ) ∫ shape a dS over the wall.
        public let volume: Volume<Space.Physical>.Value<Double>

        @inlinable
        public init(frequency: Frequency.Value<Double>, shape: [Double], volume: Volume<Space.Physical>.Value<Double>) {
            self.frequency = frequency
            self.shape = shape
            self.volume = volume
        }
    }
}

extension Acoustics.Cavity.Mode: Hashable {}
extension Acoustics.Cavity.Mode: Sendable {}
