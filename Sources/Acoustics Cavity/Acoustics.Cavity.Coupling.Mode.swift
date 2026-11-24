public import Mass

extension Acoustics.Cavity.Coupling {
    /// One mode of a cavity coupled to its walls and openings: how fast it rings, the mass that moves with it, how
    /// fast it loses its energy, and how much each wall mode and opening takes part.
    ///
    /// Struck at a point where it moves `r` per unit of its amplitude (`Coupling.reach(at:)`), it is a resonator of
    /// mass `mass / r²`, so the body's admittance there is Σ iω r² / (m (ω₀² (1 + iη) − ω²)) over its modes.
    public struct Mode {
        /// The natural frequency.
        public let frequency: Frequency.Value<Double>
        /// The modal mass at the scale where the largest amplitude of a wall mode or an opening is one.
        public let mass: Mass.Value<Double>
        /// The loss factor η: the walls' own losses, weighted by the strain energy each holds, plus the power the
        /// body radiates as a simple source into free space, ρ ω V² / (4π m c).
        public let loss: Double
        /// The volume the body pushes out at that scale: what the walls sweep plus what the openings breathe out.
        /// At low frequency the body radiates as a simple source of this volume.
        public let volume: Volume<Space.Physical>.Value<Double>
        /// Each wall's modes' amplitudes at that scale, wall by wall in the order of `Coupling.walls`.
        public let walls: [[Double]]
        /// Each opening's plug displacement out of the cavity at that scale.
        public let openings: [Double]

        @inlinable
        public init(
            frequency: Frequency.Value<Double>,
            mass: Mass.Value<Double>,
            loss: Double,
            volume: Volume<Space.Physical>.Value<Double>,
            walls: [[Double]],
            openings: [Double]
        ) {
            self.frequency = frequency
            self.mass = mass
            self.loss = loss
            self.volume = volume
            self.walls = walls
            self.openings = openings
        }
    }
}

extension Acoustics.Cavity.Coupling.Mode: Hashable {}
extension Acoustics.Cavity.Coupling.Mode: Sendable {}
