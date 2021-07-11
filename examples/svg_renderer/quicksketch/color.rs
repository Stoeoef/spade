use std::fmt::Display;

#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Copy)]
pub struct SketchColor {
    pub red: u8,
    pub green: u8,
    pub blue: u8,
}

impl Display for SketchColor {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "rgb({} {} {})", self.red, self.green, self.blue)
    }
}

impl SketchColor {
    pub fn from_rgb(red: u8, green: u8, blue: u8) -> Self {
        Self { red, green, blue }
    }

    pub const LAWN_GREEN: Self = Self {
        red: 124,
        green: 252,
        blue: 0,
    };

    pub const BLACK: Self = Self {
        red: 0,
        green: 0,
        blue: 0,
    };

    pub const DIM_GRAY: Self = Self {
        red: 105,
        green: 105,
        blue: 105,
    };
    pub const DARK_KHAKI: Self = Self {
        red: 189,
        green: 183,
        blue: 107,
    };
    pub const DODGER_BLUE: Self = Self {
        red: 30,
        green: 144,
        blue: 255,
    };
    pub const DARK_OLIVE_GREEN: Self = Self {
        red: 85,
        green: 107,
        blue: 47,
    };
    pub const CORAL: Self = Self {
        red: 255,
        green: 127,
        blue: 80,
    };
    pub const SLATE_BLUE: Self = Self {
        red: 106,
        green: 90,
        blue: 205,
    };
    pub const POWDER_BLUE: Self = Self {
        red: 176,
        green: 224,
        blue: 230,
    };
    pub const ORANGE: Self = Self {
        red: 255,
        green: 165,
        blue: 0,
    };
    pub const RED: Self = Self {
        red: 255,
        green: 0,
        blue: 0,
    };
    pub const PALE_GREEN: Self = Self {
        red: 152,
        green: 251,
        blue: 152,
    };
    pub const LIGHT_SEA_GREEN: Self = Self {
        red: 32,
        green: 178,
        blue: 170,
    };
    pub const LIGHT_GRAY: Self = Self {
        red: 211,
        green: 211,
        blue: 211,
    };
    pub const PERU: Self = Self {
        red: 205,
        green: 133,
        blue: 63,
    };
    pub const LAVENDER: Self = Self {
        red: 230,
        green: 230,
        blue: 250,
    };
    pub const INDIGO: Self = Self {
        red: 75,
        green: 0,
        blue: 130,
    };
    pub const DARK_MAGENTA: Self = Self {
        red: 139,
        green: 0,
        blue: 139,
    };
    pub const PALE_GOLDENROD: Self = Self {
        red: 238,
        green: 232,
        blue: 170,
    };
    pub const DARK_SEA_GREEN: Self = Self {
        red: 143,
        green: 188,
        blue: 143,
    };
    pub const MEDIUM_PURPLE: Self = Self {
        red: 147,
        green: 112,
        blue: 219,
    };
    pub const YELLOW_GREEN: Self = Self {
        red: 154,
        green: 205,
        blue: 50,
    };
    pub const ORCHID: Self = Self {
        red: 218,
        green: 112,
        blue: 214,
    };
    pub const PALE_TURQUOISE: Self = Self {
        red: 175,
        green: 238,
        blue: 238,
    };
    pub const MEDIUM_BLUE: Self = Self {
        red: 0,
        green: 0,
        blue: 205,
    };
    pub const FIRE_BRICK: Self = Self {
        red: 178,
        green: 34,
        blue: 34,
    };
    pub const BISQUE: Self = Self {
        red: 255,
        green: 228,
        blue: 196,
    };
    pub const PLUM: Self = Self {
        red: 221,
        green: 160,
        blue: 221,
    };
    pub const IVORY: Self = Self {
        red: 255,
        green: 255,
        blue: 240,
    };
    pub const SALMON: Self = Self {
        red: 250,
        green: 128,
        blue: 114,
    };
    pub const SLATE_GRAY: Self = Self {
        red: 112,
        green: 128,
        blue: 144,
    };
    pub const LIGHT_GREEN: Self = Self {
        red: 144,
        green: 238,
        blue: 144,
    };
    pub const VIOLET: Self = Self {
        red: 238,
        green: 130,
        blue: 238,
    };
    pub const DARK_SALMON: Self = Self {
        red: 233,
        green: 150,
        blue: 122,
    };
    pub const ALICE_BLUE: Self = Self {
        red: 240,
        green: 248,
        blue: 255,
    };
    pub const OLD_LACE: Self = Self {
        red: 253,
        green: 245,
        blue: 230,
    };
    pub const LIGHT_GOLDENROD_YELLOW: Self = Self {
        red: 250,
        green: 250,
        blue: 210,
    };
    pub const CADET_BLUE: Self = Self {
        red: 95,
        green: 158,
        blue: 160,
    };
    pub const PAPAYA_WHIP: Self = Self {
        red: 255,
        green: 239,
        blue: 213,
    };
    pub const PURPLE: Self = Self {
        red: 128,
        green: 0,
        blue: 128,
    };
    pub const LIGHT_SKY_BLUE: Self = Self {
        red: 135,
        green: 206,
        blue: 250,
    };
    pub const WHEAT: Self = Self {
        red: 245,
        green: 222,
        blue: 179,
    };
    pub const DARK_SLATE_GREY: Self = Self {
        red: 47,
        green: 79,
        blue: 79,
    };
    pub const TEAL: Self = Self {
        red: 0,
        green: 128,
        blue: 128,
    };
    pub const CORN_SILK: Self = Self {
        red: 255,
        green: 248,
        blue: 220,
    };
    pub const LIGHT_BLUE: Self = Self {
        red: 173,
        green: 216,
        blue: 230,
    };
    pub const LEMON_CHIFFON: Self = Self {
        red: 255,
        green: 250,
        blue: 205,
    };
    pub const DARK_RED: Self = Self {
        red: 139,
        green: 0,
        blue: 0,
    };
    pub const ROYAL_BLUE: Self = Self {
        red: 65,
        green: 105,
        blue: 225,
    };
    pub const MEDIUM_SEA_GREEN: Self = Self {
        red: 60,
        green: 179,
        blue: 113,
    };
    pub const ROSY_BROWN: Self = Self {
        red: 188,
        green: 143,
        blue: 143,
    };
    pub const LINEN: Self = Self {
        red: 250,
        green: 240,
        blue: 230,
    };
    pub const CORNFLOWER_BLUE: Self = Self {
        red: 100,
        green: 149,
        blue: 237,
    };
    pub const INDIAN_RED: Self = Self {
        red: 205,
        green: 92,
        blue: 92,
    };
    pub const PALE_VIOLET_RED: Self = Self {
        red: 219,
        green: 112,
        blue: 147,
    };
    pub const TURQUOISE: Self = Self {
        red: 64,
        green: 224,
        blue: 208,
    };
    pub const GOLD: Self = Self {
        red: 255,
        green: 215,
        blue: 0,
    };
    pub const AQUAMARINE: Self = Self {
        red: 127,
        green: 255,
        blue: 212,
    };
    pub const SEASHELL: Self = Self {
        red: 255,
        green: 245,
        blue: 238,
    };
    pub const DARK_TURQUOISE: Self = Self {
        red: 0,
        green: 206,
        blue: 209,
    };
    pub const HONEYDEW: Self = Self {
        red: 240,
        green: 255,
        blue: 240,
    };
    pub const TAN: Self = Self {
        red: 210,
        green: 180,
        blue: 140,
    };
    pub const ANTIQUE_WHITE: Self = Self {
        red: 250,
        green: 235,
        blue: 215,
    };
    pub const AQUA: Self = Self {
        red: 0,
        green: 255,
        blue: 255,
    };
    pub const TOMATO: Self = Self {
        red: 255,
        green: 99,
        blue: 71,
    };
    pub const LIGHT_SLATE_GRAY: Self = Self {
        red: 119,
        green: 136,
        blue: 153,
    };
    pub const GREEN: Self = Self {
        red: 0,
        green: 128,
        blue: 0,
    };
    pub const DARKBLUE: Self = Self {
        red: 0,
        green: 0,
        blue: 139,
    };
    pub const SLATE_GREY: Self = Self {
        red: 112,
        green: 128,
        blue: 144,
    };
    pub const PEACH_PUFF: Self = Self {
        red: 255,
        green: 218,
        blue: 185,
    };
    pub const DARK_SLATE_GRAY: Self = Self {
        red: 47,
        green: 79,
        blue: 79,
    };
    pub const DARK_GOLDENROD: Self = Self {
        red: 184,
        green: 134,
        blue: 11,
    };
    pub const DEEP_PINK: Self = Self {
        red: 255,
        green: 20,
        blue: 147,
    };
    pub const GREY: Self = Self {
        red: 128,
        green: 128,
        blue: 128,
    };
    pub const STEEL_BLUE: Self = Self {
        red: 70,
        green: 130,
        blue: 180,
    };
    pub const FOREST_GREEN: Self = Self {
        red: 34,
        green: 139,
        blue: 34,
    };
    pub const DARK_GRAY: Self = Self {
        red: 169,
        green: 169,
        blue: 169,
    };
    pub const LIGHT_CYAN: Self = Self {
        red: 224,
        green: 255,
        blue: 255,
    };
    pub const SILVER: Self = Self {
        red: 192,
        green: 192,
        blue: 192,
    };
    pub const BURLY_WOOD: Self = Self {
        red: 222,
        green: 184,
        blue: 135,
    };
    pub const BLUE: Self = Self {
        red: 0,
        green: 0,
        blue: 255,
    };
    pub const CYAN: Self = Self {
        red: 0,
        green: 255,
        blue: 255,
    };
    pub const SKY_BLUE: Self = Self {
        red: 135,
        green: 206,
        blue: 235,
    };
    pub const MEDIUM_AQUAMARINE: Self = Self {
        red: 102,
        green: 205,
        blue: 170,
    };
    pub const LIGHT_STEEL_BLUE: Self = Self {
        red: 176,
        green: 196,
        blue: 222,
    };
    pub const KHAKI: Self = Self {
        red: 240,
        green: 230,
        blue: 140,
    };
    pub const NAVY: Self = Self {
        red: 0,
        green: 0,
        blue: 128,
    };
    pub const PINK: Self = Self {
        red: 255,
        green: 192,
        blue: 203,
    };
    pub const BLUE_VIOLET: Self = Self {
        red: 138,
        green: 43,
        blue: 226,
    };
    pub const DARK_ORCHID: Self = Self {
        red: 153,
        green: 50,
        blue: 204,
    };
    pub const MINT_CREAM: Self = Self {
        red: 245,
        green: 255,
        blue: 250,
    };
    pub const CHOCOLATE: Self = Self {
        red: 210,
        green: 105,
        blue: 30,
    };
    pub const CHARTREUSE: Self = Self {
        red: 127,
        green: 255,
        blue: 0,
    };
    pub const LIME: Self = Self {
        red: 0,
        green: 255,
        blue: 0,
    };
    pub const MEDIUM_ORCHID: Self = Self {
        red: 186,
        green: 85,
        blue: 211,
    };
    pub const LAVENDER_BLUSH: Self = Self {
        red: 255,
        green: 240,
        blue: 245,
    };
    pub const MEDIUMSLATEBLUE: Self = Self {
        red: 123,
        green: 104,
        blue: 238,
    };
    pub const DARK_ORANGE: Self = Self {
        red: 255,
        green: 140,
        blue: 0,
    };
    pub const GHOST_WHITE: Self = Self {
        red: 248,
        green: 248,
        blue: 255,
    };
    pub const FUCHSIA: Self = Self {
        red: 255,
        green: 0,
        blue: 255,
    };
    pub const MOCCASIN: Self = Self {
        red: 255,
        green: 228,
        blue: 181,
    };
    pub const WHITE: Self = Self {
        red: 255,
        green: 255,
        blue: 255,
    };
    pub const DARK_GREY: Self = Self {
        red: 169,
        green: 169,
        blue: 169,
    };
    pub const MAROON: Self = Self {
        red: 128,
        green: 0,
        blue: 0,
    };
    pub const MIDNIGHT_BLUE: Self = Self {
        red: 25,
        green: 25,
        blue: 112,
    };
    pub const LIME_GREEN: Self = Self {
        red: 50,
        green: 205,
        blue: 50,
    };
    pub const LIGHT_CORAL: Self = Self {
        red: 240,
        green: 128,
        blue: 128,
    };
    pub const HOT_PINK: Self = Self {
        red: 255,
        green: 105,
        blue: 180,
    };
    pub const MISTY_ROSE: Self = Self {
        red: 255,
        green: 228,
        blue: 225,
    };
    pub const LIGHT_SLATE_GREY: Self = Self {
        red: 119,
        green: 136,
        blue: 153,
    };
    pub const GOLDENROD: Self = Self {
        red: 218,
        green: 165,
        blue: 32,
    };
    pub const MEDIUM_TURQUOISE: Self = Self {
        red: 72,
        green: 209,
        blue: 204,
    };
    pub const SEA_GREEN: Self = Self {
        red: 46,
        green: 139,
        blue: 87,
    };
    pub const FLORALWHITE: Self = Self {
        red: 255,
        green: 250,
        blue: 240,
    };
    pub const BLANCHEDALMOND: Self = Self {
        red: 255,
        green: 235,
        blue: 205,
    };
    pub const SPRING_GREEN: Self = Self {
        red: 0,
        green: 255,
        blue: 127,
    };
    pub const LIGHT_YELLOW: Self = Self {
        red: 255,
        green: 255,
        blue: 224,
    };
    pub const NAVAJO_WHITE: Self = Self {
        red: 255,
        green: 222,
        blue: 173,
    };
    pub const GAINSBORO: Self = Self {
        red: 220,
        green: 220,
        blue: 220,
    };
    pub const GREEN_YELLOW: Self = Self {
        red: 173,
        green: 255,
        blue: 47,
    };
    pub const DEEP_SKY_BLUE: Self = Self {
        red: 0,
        green: 191,
        blue: 255,
    };
    pub const SANDY_BROWN: Self = Self {
        red: 244,
        green: 164,
        blue: 96,
    };
    pub const AZURE: Self = Self {
        red: 240,
        green: 255,
        blue: 255,
    };
    pub const BROWN: Self = Self {
        red: 165,
        green: 42,
        blue: 42,
    };
    pub const MAGENTA: Self = Self {
        red: 255,
        green: 0,
        blue: 255,
    };
    pub const DIM_GREY: Self = Self {
        red: 105,
        green: 105,
        blue: 105,
    };
    pub const MEDIUM_VIOLET_RED: Self = Self {
        red: 199,
        green: 21,
        blue: 133,
    };
    pub const SNOW: Self = Self {
        red: 255,
        green: 250,
        blue: 250,
    };
    pub const YELLOW: Self = Self {
        red: 255,
        green: 255,
        blue: 0,
    };
    pub const GRAY: Self = Self {
        red: 128,
        green: 128,
        blue: 128,
    };
    pub const ORANGE_RED: Self = Self {
        red: 255,
        green: 69,
        blue: 0,
    };
    pub const CRIMSON: Self = Self {
        red: 220,
        green: 20,
        blue: 60,
    };
    pub const DARK_GREEN: Self = Self {
        red: 0,
        green: 100,
        blue: 0,
    };
    pub const LIGHT_PINK: Self = Self {
        red: 255,
        green: 182,
        blue: 193,
    };
    pub const MEDIUM_SPRING_GREEN: Self = Self {
        red: 0,
        green: 250,
        blue: 154,
    };
    pub const THISTLE: Self = Self {
        red: 216,
        green: 191,
        blue: 216,
    };
    pub const SIENNA: Self = Self {
        red: 160,
        green: 82,
        blue: 45,
    };
    pub const LIGHT_GREY: Self = Self {
        red: 211,
        green: 211,
        blue: 211,
    };
    pub const DARK_SLATE_BLUE: Self = Self {
        red: 72,
        green: 61,
        blue: 139,
    };
    pub const LIGHT_SALMON: Self = Self {
        red: 255,
        green: 160,
        blue: 122,
    };
    pub const DARK_VIOLET: Self = Self {
        red: 148,
        green: 0,
        blue: 211,
    };
    pub const SADDLE_BROWN: Self = Self {
        red: 139,
        green: 69,
        blue: 19,
    };
    pub const DARK_CYAN: Self = Self {
        red: 0,
        green: 139,
        blue: 139,
    };
    pub const OLIVE: Self = Self {
        red: 128,
        green: 128,
        blue: 0,
    };
    pub const WHITE_SMOKE: Self = Self {
        red: 245,
        green: 245,
        blue: 245,
    };
    pub const BEIGE: Self = Self {
        red: 245,
        green: 245,
        blue: 220,
    };
    pub const OLIVE_DRAB: Self = Self {
        red: 107,
        green: 142,
        blue: 35,
    };
}
