extension Acoustics.Cavity {
    /// A hole in a closed cavity's wall, such as a guitar's sound hole or a lute's rose: the plug of air in it moves
    /// as one piston against the pressure of the whole cavity, as in Helmholtz's resonator.
    ///
    /// The plug's mass is ρ S L, its length L including the end corrections on both sides, 0.85 a each for a
    /// circular hole of radius a in a wide flat wall (Ingard 1953). It is driven by the cavity's uniform pressure
    /// only, so the cavity must be closed and its zero mode among the modes coupled.
    public struct Opening {
        /// The area S of the hole.
        public let area: Area<Space.Physical>.Value<Double>
        /// The effective length L of the plug: the wall's thickness plus both end corrections.
        public let length: Measure<1, Space.Physical>.Value<Double>

        @inlinable
        public init(area: Area<Space.Physical>.Value<Double>, length: Measure<1, Space.Physical>.Value<Double>) {
            self.area = area
            self.length = length
        }
    }
}

extension Acoustics.Cavity.Opening {
    /// A circular hole of `radius` in a flat wall of `thickness`: its area π a², and its length the thickness plus
    /// 0.85 a on each side (Ingard 1953, J. Acoust. Soc. Am. 25, 1037, for a flanged end).
    @inlinable
    public static func circle(
        radius: Measure<1, Space.Physical>.Value<Double>,
        thickness: Measure<1, Space.Physical>.Value<Double>
    ) -> Self {
        let a = radius.underlying
        return Self(
            area: Area<Space.Physical>.amount(Double.pi * a * a),
            length: Measure<1, Space.Physical>.amount(thickness.underlying + 2 * 0.85 * a)
        )
    }
}

extension Acoustics.Cavity.Opening: Hashable {}
extension Acoustics.Cavity.Opening: Sendable {}
