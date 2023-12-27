use clap::Parser;
use num_complex::ComplexFloat;
use num_complex::{Complex32, Complex64};
use reed_solomon::Decoder;
use rustfft::FftPlanner;
use std::collections::vec_deque::VecDeque;
use std::fs;
use std::io;
use std::io::Read;
use std::io::Write;

const TS_SYNC_BYTE: u8 = 0x47;
const TS_SIZE: usize = 188;
const TS_PARITY_SIZE: usize = 16;
const TSP_SIZE: usize = TS_SIZE + TS_PARITY_SIZE;
const BYTE_INTERLEAVER_SIZE: usize = 12;

const AC_CARRIER: &[&[usize]] = &[
    &[10, 28, 161, 191, 277, 316, 335, 425],
    &[20, 40, 182, 208, 251, 295, 400, 421],
    &[4, 89, 148, 197, 224, 280, 331, 413],
    &[98, 101, 118, 136, 269, 299, 385, 424],
    &[11, 101, 128, 148, 290, 316, 359, 403],
    &[76, 97, 112, 197, 256, 305, 332, 388],
    &[7, 89, 206, 209, 226, 244, 377, 407],
    &[61, 100, 119, 209, 236, 256, 398, 424],
    &[35, 79, 184, 205, 220, 305, 364, 413],
    &[8, 64, 115, 197, 314, 317, 334, 352],
    &[53, 83, 169, 208, 227, 317, 344, 364],
    &[74, 100, 143, 187, 292, 313, 328, 413],
    &[40, 89, 116, 172, 223, 305, 422, 425],
];

const TMCC_CARRIER: &[&[usize]] = &[
    &[70, 133, 233, 410],
    &[44, 155, 265, 355],
    &[83, 169, 301, 425],
    &[23, 178, 241, 341],
    &[86, 152, 263, 373],
    &[31, 191, 277, 409],
    &[101, 131, 286, 349],
    &[17, 194, 260, 371],
    &[49, 139, 299, 385],
    &[85, 209, 239, 394],
    &[25, 125, 302, 368],
    &[47, 157, 247, 407],
    &[61, 193, 317, 347],
];

const MODE1_CARRIER_RANDOMIZE: &[usize] = &[
    80, 93, 63, 92, 94, 55, 17, 81, 6, 51, 9, 85, 89, 65, 52, 15, 73, 66, 46, 71, 12, 70, 18, 13,
    95, 34, 1, 38, 78, 59, 91, 64, 0, 28, 11, 4, 45, 35, 16, 7, 48, 22, 23, 77, 56, 19, 8, 36, 39,
    61, 21, 3, 26, 69, 67, 20, 74, 86, 72, 25, 31, 5, 49, 42, 54, 87, 43, 60, 29, 2, 76, 84, 83,
    40, 14, 79, 27, 57, 44, 37, 30, 68, 47, 88, 75, 41, 90, 10, 33, 32, 62, 50, 58, 82, 53, 24,
];

const MODE2_CARRIER_RANDOMIZE: &[usize] = &[
    98, 35, 67, 116, 135, 17, 5, 93, 73, 168, 54, 143, 43, 74, 165, 48, 37, 69, 154, 150, 107, 76,
    176, 79, 175, 36, 28, 78, 47, 128, 94, 163, 184, 72, 142, 2, 86, 14, 130, 151, 114, 68, 46,
    183, 122, 112, 180, 42, 105, 97, 33, 134, 177, 84, 170, 45, 187, 38, 167, 10, 189, 51, 117,
    156, 161, 25, 89, 125, 139, 24, 19, 57, 71, 39, 77, 191, 88, 85, 0, 162, 181, 113, 140, 61, 75,
    82, 101, 174, 118, 20, 136, 3, 121, 190, 120, 92, 160, 52, 153, 127, 65, 60, 133, 147, 131, 87,
    22, 58, 100, 111, 141, 83, 49, 132, 12, 155, 146, 102, 164, 66, 1, 62, 178, 15, 182, 96, 80,
    119, 23, 6, 166, 56, 99, 123, 138, 137, 21, 145, 185, 18, 70, 129, 95, 90, 149, 109, 124, 50,
    11, 152, 4, 31, 172, 40, 13, 32, 55, 159, 41, 8, 7, 144, 16, 26, 173, 81, 44, 103, 64, 9, 30,
    157, 126, 179, 148, 63, 188, 171, 106, 104, 158, 115, 34, 186, 29, 108, 53, 91, 169, 110, 27,
    59,
];

const MODE3_CARRIER_RANDOMIZE: &[usize] = &[
    62, 13, 371, 11, 285, 336, 365, 220, 226, 92, 56, 46, 120, 175, 298, 352, 172, 235, 53, 164,
    368, 187, 125, 82, 5, 45, 173, 258, 135, 182, 141, 273, 126, 264, 286, 88, 233, 61, 249, 367,
    310, 179, 155, 57, 123, 208, 14, 227, 100, 311, 205, 79, 184, 185, 328, 77, 115, 277, 112, 20,
    199, 178, 143, 152, 215, 204, 139, 234, 358, 192, 309, 183, 81, 129, 256, 314, 101, 43, 97,
    324, 142, 157, 90, 214, 102, 29, 303, 363, 261, 31, 22, 52, 305, 301, 293, 177, 116, 296, 85,
    196, 191, 114, 58, 198, 16, 167, 145, 119, 245, 113, 295, 193, 232, 17, 108, 283, 246, 64, 237,
    189, 128, 373, 302, 320, 239, 335, 356, 39, 347, 351, 73, 158, 276, 243, 99, 38, 287, 3, 330,
    153, 315, 117, 289, 213, 210, 149, 383, 337, 339, 151, 241, 321, 217, 30, 334, 161, 322, 49,
    176, 359, 12, 346, 60, 28, 229, 265, 288, 225, 382, 59, 181, 170, 319, 341, 86, 251, 133, 344,
    361, 109, 44, 369, 268, 257, 323, 55, 317, 381, 121, 360, 260, 275, 190, 19, 63, 18, 248, 9,
    240, 211, 150, 230, 332, 231, 71, 255, 350, 355, 83, 87, 154, 218, 138, 269, 348, 130, 160,
    278, 377, 216, 236, 308, 223, 254, 25, 98, 300, 201, 137, 219, 36, 325, 124, 66, 353, 169, 21,
    35, 107, 50, 106, 333, 326, 262, 252, 271, 263, 372, 136, 0, 366, 206, 159, 122, 188, 6, 284,
    96, 26, 200, 197, 186, 345, 340, 349, 103, 84, 228, 212, 2, 67, 318, 1, 74, 342, 166, 194, 33,
    68, 267, 111, 118, 140, 195, 105, 202, 291, 259, 23, 171, 65, 281, 24, 165, 8, 94, 222, 331,
    34, 238, 364, 376, 266, 89, 80, 253, 163, 280, 247, 4, 362, 379, 290, 279, 54, 78, 180, 72,
    316, 282, 131, 207, 343, 370, 306, 221, 132, 7, 148, 299, 168, 224, 48, 47, 357, 313, 75, 104,
    70, 147, 40, 110, 374, 69, 146, 37, 375, 354, 174, 41, 32, 304, 307, 312, 15, 272, 134, 242,
    203, 209, 380, 162, 297, 327, 10, 93, 42, 250, 156, 338, 292, 144, 378, 294, 329, 127, 270, 76,
    95, 91, 244, 274, 27, 51,
];

fn calc_mmse(signal: &[Complex32], symbol_len: usize, result: &mut [f64]) {
    for i in 0..signal.len() - symbol_len {
        let c1 = Complex64::new(signal[i].re as f64, signal[i].im as f64).norm_sqr();
        let c2 = Complex64::new(
            signal[i + symbol_len].re as f64,
            signal[i + symbol_len].im as f64,
        )
        .norm_sqr();
        let c3 = Complex64::new(signal[i].re as f64, signal[i].im as f64)
            * Complex64::new(
                signal[i + symbol_len].re as f64,
                signal[i + symbol_len].im as f64,
            )
            .conj();
        result[i] = c1 + c2 - 2f64 * c3.norm();
    }
}

fn detect_symbol_timing_metrics(
    guard_interval_len: usize,
    search_range: usize,
    metrics: &[f64],
) -> usize {
    let mut c: f64 = 0.0;
    for j in 0..guard_interval_len {
        c += metrics[j];
    }
    let mut min: f64 = c;
    let mut min_i = 0usize;
    for i in 1..search_range {
        c -= metrics[i - 1];
        c += metrics[i + guard_interval_len - 1];
        if min > c {
            min = c;
            min_i = i;
        }
    }
    return min_i;
}

fn estimate_carrier_offset(signal: &[Complex32], symbol_len: usize, guard_interval_len: usize) -> f32 {
    let mut r = Complex32::default();
    for i in 0..guard_interval_len {
        r +=signal[i].conj() * signal[i + symbol_len];
    }
    return r.arg();
}

fn read_iq(
    reader: &mut Box<dyn io::Read>,
    buffer: &mut [u8],
    iq_buffer: &mut [Complex32],
    bit: Bit,
) -> io::Result<()> {
    match bit {
        Bit::Bit16 => {
            let iq_size = 4;
            reader.read_exact(&mut buffer[..iq_size * iq_buffer.len()])?;
            for i in 0..iq_buffer.len() {
                let i_component: i16 = i16::from_le_bytes(
                    buffer[i * iq_size..i * iq_size + iq_size / 2]
                        .try_into()
                        .unwrap(),
                );
                let q_component: i16 = i16::from_le_bytes(
                    buffer[i * iq_size + iq_size / 2..i * iq_size + iq_size]
                        .try_into()
                        .unwrap(),
                );
                iq_buffer[i] = Complex32::new(i_component as f32, q_component as f32);
            }
        }
        Bit::Bit8 => {
            let iq_size = 2;
            reader.read_exact(&mut buffer[..iq_size * iq_buffer.len()])?;
            for i in 0..iq_buffer.len() {
                let i_component: i8 = i8::from_le_bytes(
                    buffer[i * iq_size..i * iq_size + iq_size / 2]
                        .try_into()
                        .unwrap(),
                );
                let q_component: i8 = i8::from_le_bytes(
                    buffer[i * iq_size + iq_size / 2..i * iq_size + iq_size]
                        .try_into()
                        .unwrap(),
                );
                iq_buffer[i] = Complex32::new(i_component as f32, q_component as f32);
            }
        }
    }
    return io::Result::Ok(());
}

fn estimate_integer_carrier_offset(
    f: &[Complex32],
    min: usize,
    max: usize,
    cor_carriers: &[usize],
) -> usize {
    let mut max_g = f32::MIN;
    let mut max_i = usize::MAX;
    for i in min..=max {
        let mut g = 0f32;
        for c in cor_carriers {
            g += f[i + c].norm_sqr();
        }
        if g > max_g {
            max_g = g;
            max_i = i;
        }
    }
    return max_i;
}

fn estimate_scattered_pilot(f: &[Complex32]) -> usize {
    let mut max_g = f32::MIN;
    let mut max_i = usize::MAX;
    for pos in 0..4 {
        let mut g = 0f32;
        for c in 0..f.len() / 12 {
            let s = f[pos * 3 + c * 12].norm_sqr();
            assert_ne!(s, 0.0f32);
            g += s;
        }
        if g > max_g {
            max_g = g;
            max_i = pos * 3;
        }
    }
    return max_i;
}

const PILOT_PRBS_INITIAL_STATE: u32 = 0b11111111111;

fn pilot_prbs(d: u32) -> u32 {
    let d11 = d & 1;
    let d9 = if (d & 0b100) == 0b100 { 1 } else { 0 };
    return (d >> 1) | ((d11 ^ d9) << 10);
}

const BYTE_PRBS_INITIAL_STATE: u16 = 0b100101010000000;

fn byte_prbs(mut d: u16) -> u16 {
    let d15 = d & 1;
    let d14 = (d & 2) >> 1;
    d >>= 1;
    if (d15 ^ d14) != 0 {
        d |= 0b100000000000000;
    }
    return d;
}

fn bpsk_symbol(b: u32) -> Complex32 {
    if b == 0 {
        return Complex32::new(4.0 / 3.0, 0.0);
    } else {
        return Complex32::new(-4.0 / 3.0, 0.0);
    }
}

fn demap_bpsk(symbol: Complex32) -> u8 {
    if symbol.re >= 0.0 {
        return 0;
    } else {
        return 1;
    }
}

fn demodulate_dbpsk(dbpsk: &[u8]) -> Vec<u8> {
    let mut i = dbpsk[0];
    let mut data = Vec::<u8>::with_capacity((dbpsk.len() - 1 + 7) / 8);
    data.resize((dbpsk.len() - 1 + 7) / 8, 0);
    for (pos, b) in dbpsk.iter().skip(1).enumerate() {
        if (i ^ b) & 1 == 1 {
            data[pos / 8] |= 1u8 << (pos % 8);
        }
        i = *b;
    }
    return data;
}

fn equalize(segment_carriers: &mut [Complex32], equalizer: &[Complex32]) {
    for (s, e) in segment_carriers.iter_mut().zip(equalizer) {
        *s /= e;
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Copy)]
enum CodingRate {
    Rate1_2,
    Rate2_3,
    Rate3_4,
    Rate5_6,
    Rate7_8,
}

fn depuncture(data: &[f32], rate: CodingRate) -> Vec<(f32, f32)> {
    match rate {
        CodingRate::Rate1_2 => {
            let mut output = Vec::<(f32, f32)>::with_capacity(data.len() * 1 / 2);
            for i in 0..data.len() / 2 {
                output.push((data[i * 2], data[i * 2 + 1]));
            }
            return output;
        }
        CodingRate::Rate2_3 => {
            let mut output = Vec::<(f32, f32)>::with_capacity(data.len() * 2 / 3);
            for i in 0..data.len() / 3 {
                output.push((data[i * 3], data[i * 3 + 1]));
                output.push((0.0, data[i * 3 + 2]));
            }
            return output;
        }
        CodingRate::Rate3_4 => {
            let mut output = Vec::<(f32, f32)>::with_capacity(data.len() * 3 / 4);
            for i in 0..data.len() / 4 {
                output.push((data[i * 4], data[i * 4 + 1]));
                output.push((0.0, data[i * 4 + 2]));
                output.push((data[i * 4 + 3], 0.0));
            }
            return output;
        }
        CodingRate::Rate5_6 => {
            let mut output = Vec::<(f32, f32)>::with_capacity(data.len() * 5 / 6);
            for i in 0..data.len() / 6 {
                output.push((data[i * 6], data[i * 6 + 1]));
                output.push((0.0, data[i * 6 + 2]));
                output.push((data[i * 6 + 3], 0.0));
                output.push((0.0, data[i * 6 + 4]));
                output.push((data[i * 6 + 5], 0.0));
            }
            return output;
        }
        CodingRate::Rate7_8 => {
            let mut output = Vec::<(f32, f32)>::with_capacity(data.len() * 7 / 8);
            for i in 0..data.len() / 8 {
                output.push((data[i * 8], data[i * 8 + 1]));
                output.push((0.0, data[i * 8 + 2]));
                output.push((0.0, data[i * 8 + 3]));
                output.push((0.0, data[i * 8 + 4]));
                output.push((data[i * 8 + 5], 0.0));
                output.push((0.0, data[i * 8 + 6]));
                output.push((data[i * 8 + 7], 0.0));
            }
            return output;
        }
    }
}

pub fn viterbi(conv: &[(f32, f32)], initial_state: u8) -> Vec<u8> {
    const G1: u8 = 0o171;
    const G2: u8 = 0o133;
    let max_path = 1 << 7;
    let path_length = 31;
    let mut path_list = Vec::<Path>::with_capacity(max_path);
    let mut next_path_list = Vec::<Path>::with_capacity(max_path);
    path_list.resize(
        max_path,
        Path {
            d: 0,
            distance: f32::INFINITY,
            decoded: 0,
        },
    );
    next_path_list.resize(
        max_path,
        Path {
            d: 0,
            distance: f32::INFINITY,
            decoded: 0,
        },
    );
    let mut decoded_data = Vec::<u8>::with_capacity((conv.len() + 7) / 8);
    decoded_data.resize((conv.len() + 7) / 8, 0);
    #[derive(Debug, Clone)]
    struct Path {
        d: u8,
        distance: f32,
        decoded: u32,
    }

    path_list[initial_state as usize] = Path {
        d: initial_state,
        distance: 0.0,
        decoded: 0,
    };
    let mut decoded_pos = 0usize;
    let mut min_distance = 0f32;
    for (pos, (x, y)) in conv.iter().enumerate() {
        next_path_list.fill(Path {
            d: 0,
            distance: f32::INFINITY,
            decoded: 0,
        });
        for path in &path_list {
            if path.distance.is_infinite() {
                continue;
            }
            // if 0
            let d_0 = path.d >> 1;
            let x_output_0 = if (d_0 & G1).count_ones() % 2 == 1 {
                -1.0f32
            } else {
                1.0f32
            };
            let y_output_0 = if (d_0 & G2).count_ones() % 2 == 1 {
                -1.0f32
            } else {
                1.0f32
            };
            let distance_0 =
                path.distance + (x_output_0 - x).powi(2) + (y_output_0 - y).powi(2) - min_distance;
            if next_path_list[d_0 as usize].distance > distance_0 {
                next_path_list[d_0 as usize] = Path {
                    d: d_0,
                    distance: distance_0,
                    decoded: path.decoded << 1,
                };
            }
            // if 1
            let d_1 = (path.d >> 1) | 0b1000000;
            let x_output_1 = if (d_1 & G1).count_ones() % 2 == 1 {
                -1.0f32
            } else {
                1.0f32
            };
            let y_output_1 = if (d_1 & G2).count_ones() % 2 == 1 {
                -1.0f32
            } else {
                1.0f32
            };
            let distance_1 =
                path.distance + (x_output_1 - x).powi(2) + (y_output_1 - y).powi(2) - min_distance;
            if next_path_list[d_1 as usize].distance > distance_1 {
                next_path_list[d_1 as usize] = Path {
                    d: d_1,
                    distance: distance_1,
                    decoded: 1 | (path.decoded << 1),
                };
            }
        }
        std::mem::swap(&mut path_list, &mut next_path_list);
        if pos >= path_length {
            let min = path_list
                .iter()
                .min_by(|a, b| a.distance.partial_cmp(&b.distance).unwrap())
                .unwrap();
            min_distance = min.distance;
            let decoded = min.decoded;
            if (decoded & (1 << path_length)) != 0 {
                decoded_data[decoded_pos / 8] |= 1 << (7 - (decoded_pos % 8));
            }
            decoded_pos += 1;
        }
    }
    let mut decoded = path_list
        .iter()
        .min_by(|a, b| a.distance.partial_cmp(&b.distance).unwrap())
        .unwrap()
        .decoded;
    for _ in decoded_pos..conv.len() {
        decoded <<= 1;
        if (decoded & (1 << path_length)) != 0 {
            decoded_data[decoded_pos / 8] |= 1 << (7 - (decoded_pos % 8));
        }
        decoded_pos += 1;
    }
    return decoded_data;
}

const INVALID_COMPLEX32: Complex32 = Complex32::new(f32::INFINITY, f32::INFINITY);

#[derive(Debug, Clone, PartialEq, Eq, Copy)]
pub enum Mode {
    Mode1,
    Mode2,
    Mode3,
}

pub struct TimeDeinterleaver {
    buffers: Vec<VecDeque<Complex32>>,
}

impl TimeDeinterleaver {
    pub fn new(mode: Mode, time_interleaving_length: u8) -> TimeDeinterleaver {
        let i = match (time_interleaving_length, mode) {
            (0, _) => 0,
            (n, Mode::Mode1) => (1 << (n - 1)) * 4,
            (n, Mode::Mode2) => (1 << (n - 1)) * 2,
            (n, Mode::Mode3) => (1 << (n - 1)) * 1,
        };
        let carriers = match mode {
            Mode::Mode1 => MODE1_CARRIER_RANDOMIZE.len(),
            Mode::Mode2 => MODE2_CARRIER_RANDOMIZE.len(),
            Mode::Mode3 => MODE3_CARRIER_RANDOMIZE.len(),
        };
        let mut buffers = Vec::<VecDeque<Complex32>>::with_capacity(carriers);
        let max_delay = (0..carriers).map(|j| i * ((j * 5) % 96)).max().unwrap();
        for j in 0..carriers {
            let size = max_delay - (i * ((j * 5) % 96));
            let mut symbol_buffer = VecDeque::<Complex32>::with_capacity(size + 1);
            symbol_buffer.resize(size, INVALID_COMPLEX32);
            buffers.push(symbol_buffer);
        }
        return TimeDeinterleaver { buffers };
    }
    pub fn push(&mut self, carrier: usize, symbol: Complex32) {
        self.buffers[carrier].push_back(symbol);
    }
    pub fn pop(&mut self, carrier: usize) -> Complex32 {
        return self.buffers[carrier].pop_front().unwrap();
    }
}

pub struct QPSKDeinterleaver {
    b0: VecDeque<f32>,
    b1: VecDeque<f32>,
}

impl QPSKDeinterleaver {
    fn new() -> QPSKDeinterleaver {
        let mut b0 = VecDeque::<f32>::with_capacity(120 + 1);
        b0.resize(120, f32::INFINITY);
        let mut b1 = VecDeque::<f32>::with_capacity(0 + 1);
        b1.resize(0, f32::INFINITY);
        return QPSKDeinterleaver { b0, b1 };
    }
    pub fn push(&mut self, symbol: Complex32) {
        self.b0.push_back(symbol.re);
        self.b1.push_back(symbol.im);
    }
    pub fn pop(&mut self) -> (f32, f32) {
        return (self.b0.pop_front().unwrap(), self.b1.pop_front().unwrap());
    }
}

pub struct ByteDeinterleaver {
    buffers: Vec<(usize, VecDeque<u8>)>,
}

impl ByteDeinterleaver {
    pub fn new() -> ByteDeinterleaver {
        let mut buffers = Vec::<(usize, VecDeque<u8>)>::with_capacity(BYTE_INTERLEAVER_SIZE);
        for i in 0..BYTE_INTERLEAVER_SIZE {
            let size = 17 * (BYTE_INTERLEAVER_SIZE - 1) - 17 * i;
            let mut d = VecDeque::with_capacity(size + 1);
            d.resize(size, 0);
            buffers.push((size + 1, d));
        }
        return ByteDeinterleaver { buffers };
    }
    pub fn push(&mut self, i: usize, byte: u8) {
        self.buffers[i % BYTE_INTERLEAVER_SIZE].1.push_back(byte);
    }
    pub fn pop(&mut self, i: usize) -> Option<u8> {
        let (size, buf) = &mut self.buffers[i % BYTE_INTERLEAVER_SIZE];
        let b = buf.pop_front().unwrap();
        if *size != 1 {
            *size -= 1;
            return None;
        }
        return Some(b);
    }
}

fn convolve(data: &[u8]) -> Vec<(f32, f32)> {
    let mut v = Vec::with_capacity(data.len());
    let mut d = 0u8;
    const G1: u8 = 0o171;
    const G2: u8 = 0o133;
    for b in data {
        for i in (0..8).rev() {
            if (b & (1 << i)) != 0 {
                d |= 0b1000000;
            }
            let x = if (d & G1).count_ones() % 2 == 1 {
                -1.0f32
            } else {
                1.0f32
            };
            let y = if (d & G2).count_ones() % 2 == 1 {
                -1.0f32
            } else {
                1.0f32
            };
            v.push((x, y));
            d >>= 1;
        }
    }
    return v;
}

fn estimate_cnr(decoded: &[(f32, f32)], expected: &[(f32, f32)]) -> f64 {
    let mut noise_power = 0.0f64;
    let mut signal_power = 0.0f64;
    for ((dx, dy), (ex, ey)) in decoded.iter().zip(expected.iter()) {
        if *dx != 0.0f32 {
            noise_power += (dx - ex).powi(2) as f64;
            signal_power += 1.0f64;
        }
        if *dy != 0.0f32 {
            noise_power += (dy - ey).powi(2) as f64;
            signal_power += 1.0f64;
        }
    }
    return 20.0f64 * (signal_power.sqrt() / noise_power.sqrt()).log10();
}

#[derive(Debug, Clone, clap::ValueEnum, Copy)]
enum Bit {
    #[clap(name = "16")]
    Bit16,
    #[clap(name = "8")]
    Bit8,
}
#[derive(Parser, Debug)]
#[command(author, version, about)]
struct Args {
    #[arg(short, long)]
    input: Option<String>,
    #[arg(short, long)]
    output: String,
    #[arg(short, long)]
    verbose: bool,
    #[arg(value_enum, long, default_value_t = Bit::Bit16)]
    bit: Bit,
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = Args::parse();
    let mut reader: Box<dyn io::Read> = if let Some(input) = &args.input {
        Box::new(io::BufReader::new(
            fs::OpenOptions::new().read(true).open(input)?,
        ))
    } else {
        Box::new(io::BufReader::new(std::io::stdin()))
    };
    let out_file = fs::OpenOptions::new()
        .truncate(true)
        .create(true)
        .write(true)
        .open(&args.output)?;
    let mut writer = io::BufWriter::with_capacity(204 * 40, out_file);
    let symbol_len = 8192;
    let guard_interval_len = symbol_len / 8;
    let buf_size: usize = (guard_interval_len + symbol_len) * 2;
    let iq_size: usize = 4;
    let mut buffer = Vec::<u8>::with_capacity(buf_size * iq_size);
    let mut iq_buffer = Vec::<Complex32>::with_capacity(buf_size);
    buffer.resize(buf_size * iq_size, 0);
    iq_buffer.resize(buf_size, Complex32::default());
    reader.read_exact(&mut buffer)?;
    read_iq(&mut reader, &mut buffer, &mut iq_buffer, args.bit)?;
    let mut planner = FftPlanner::<f32>::new();
    let fft = planner.plan_fft_forward(symbol_len);
    let mut carriers = Vec::<Complex32>::with_capacity(symbol_len);
    carriers.resize(symbol_len, Complex32::default());
    let symbols_per_segment = 432usize;
    let mut ccc = Vec::<usize>::new();
    for (i, ac) in AC_CARRIER.iter().enumerate() {
        for c in *ac {
            ccc.push(i * symbols_per_segment + c);
        }
    }
    for (i, ac) in TMCC_CARRIER.iter().enumerate() {
        for c in *ac {
            ccc.push(i * symbols_per_segment + c);
        }
    }
    let mut prbs_table = Vec::<u32>::new();
    let mut prbs_state = PILOT_PRBS_INITIAL_STATE;
    for _ in 0..symbols_per_segment * 13 {
        prbs_table.push(prbs_state & 1);
        prbs_state = pilot_prbs(prbs_state);
    }
    let mut equalizer = Vec::<Complex32>::with_capacity(symbols_per_segment);
    equalizer.resize(symbols_per_segment, Complex32::default());
    let mut tmcc = Vec::<Vec<u8>>::new();
    let frame_size = 204;
    tmcc.resize_with(TMCC_CARRIER[0].len(), || {
        Vec::<u8>::with_capacity(frame_size)
    });
    let mut mode3_carrier_derandomize = Vec::<usize>::with_capacity(MODE3_CARRIER_RANDOMIZE.len());
    mode3_carrier_derandomize.resize(MODE3_CARRIER_RANDOMIZE.len(), 0);
    for (i, carrier) in MODE3_CARRIER_RANDOMIZE.iter().enumerate() {
        mode3_carrier_derandomize[*carrier] = i;
    }
    let mode = Mode::Mode3;
    let mut layer_a = LayerParameter {
        carrier_modulation: TMCCData::Data(CarrierModulation::QPSK),
        coding_rate: TMCCData::Data(CodingRate::Rate2_3),
        number_of_segments: TMCCData::Data(1),
        time_interleaving_length: TMCCData::Data(0b11),
    };
    let mut time_deinterleaver =
        TimeDeinterleaver::new(mode, layer_a.time_interleaving_length.unwrap());
    let mut frame_symbols = Vec::<f32>::new();
    let mut bit_deinterleaver = QPSKDeinterleaver::new();
    let mut byte_deinterleaver = ByteDeinterleaver::new();
    let mut prev_sp_offset = 100;
    let mut frame_number = 0;
    let search_range = symbol_len;
    let mut metrics = Vec::<f64>::with_capacity(symbol_len * 2);
    metrics.resize(symbol_len * 2, 0.0);
    let mut prev_cfo = usize::MAX;
    let cnr_interval = 10000;
    let mut cnr_counter = 0;
    let mut offset = 0usize;
    loop {
        calc_mmse(
            &iq_buffer[0..search_range * 2 + guard_interval_len],
            search_range,
            &mut metrics,
        );
        let m: usize = detect_symbol_timing_metrics(guard_interval_len, search_range, &metrics);
        if args.verbose {
            eprintln!("valid symbol offset: {m}");
        }
        let delta = estimate_carrier_offset(&iq_buffer[m..], symbol_len, guard_interval_len);
        let t_s = 1.0f32 / 8126984.0f32;
        let delta = delta / (2.0f32 * std::f32::consts::PI * symbol_len as f32 * t_s);
        if args.verbose {
            eprintln!("freq offset: {delta}");
        }
        let symbol_begin = m + guard_interval_len;
        let symbol_end = m + guard_interval_len + symbol_len;
        let next_off = m + symbol_len;
        carriers.copy_from_slice(&iq_buffer[symbol_begin..symbol_end]);
        let mut dc_offset = Complex32::default();
        for c in &carriers {
            dc_offset += c;
        }
        dc_offset /= carriers.len() as f32;
        for (i, c) in carriers.iter_mut().enumerate() {
            let t = (offset + i + symbol_begin) as f32 * t_s;
            *c = (*c - dc_offset) * Complex32::new(0.0, -2.0f32 * std::f32::consts::PI * delta * t).exp();
        }
        fft.process(&mut carriers);
        carriers[0] = Complex32::default(); // remove DC offset
        let (positive_freq, negative_freq) = carriers.split_at_mut(symbol_len / 2);
        positive_freq.swap_with_slice(negative_freq);
        let estimate_range = 80;
        let carrier_offset = symbol_len / 2 - symbols_per_segment / 2 - symbols_per_segment * 6;
        let off = estimate_integer_carrier_offset(
            &carriers,
            carrier_offset - estimate_range,
            carrier_offset + estimate_range,
            &ccc,
        );
        if prev_cfo != off {
            if prev_cfo != usize::MAX {
                eprintln!(
                    "CFO changed {} -> {}",
                    prev_cfo as isize - carrier_offset as isize,
                    off as isize - carrier_offset as isize
                );
            } else {
                eprintln!("CFO: {}", off as isize - carrier_offset as isize);
            }
        }
        prev_cfo = off;
        let first_segment = off;
        let segment = 6; // 1seg
        let segment_carriers = &mut carriers[first_segment + symbols_per_segment * segment
            ..first_segment + symbols_per_segment * (segment + 1)];
        let segment_prbs =
            &prbs_table[symbols_per_segment * segment..symbols_per_segment * (segment + 1)];
        let sp_offset = estimate_scattered_pilot(&segment_carriers);
        if prev_sp_offset != 100 && (prev_sp_offset + 3) % 12 != sp_offset {
            eprintln!("SP offset {} != {sp_offset}", (prev_sp_offset + 3) % 12);
        }
        prev_sp_offset = sp_offset;
        for i in 0..segment_carriers.len() / 12 {
            equalizer[i * 12 + sp_offset] = segment_carriers[i * 12 + sp_offset]
                / bpsk_symbol(segment_prbs[i * 12 + sp_offset]);
        }
        for i in 0..segment_carriers.len() / 12 - 1 {
            let index = i * 12 + sp_offset;
            let next_index = index + 12;
            let e1 = equalizer[index];
            let e2 = equalizer[next_index];
            for j in 1..12 {
                let a = j as f32 / 12f32;
                equalizer[index + j] = e1 * (1f32 - a) + e2 * a;
            }
        }
        for j in 0..sp_offset {
            equalizer[j] = equalizer[sp_offset];
        }
        let last_sp = segment_carriers.len() / 12 * 12 - 12 + sp_offset;
        for j in last_sp + 1..segment_carriers.len() {
            equalizer[j] = equalizer[last_sp];
        }
        equalize(segment_carriers, &equalizer);
        for (carrier, i) in (0..symbols_per_segment)
            .filter(|i| {
                i % 12 != sp_offset
                    && !TMCC_CARRIER[segment].contains(i)
                    && !AC_CARRIER[segment].contains(i)
            })
            .enumerate()
        {
            time_deinterleaver.push(mode3_carrier_derandomize[carrier], segment_carriers[i]);
        }
        for carrier in 0..MODE3_CARRIER_RANDOMIZE.len() {
            let s = time_deinterleaver.pop(carrier);
            bit_deinterleaver.push(s * 2.0f32.sqrt());
            let (b0, b1) = bit_deinterleaver.pop();
            frame_symbols.push(b0);
            frame_symbols.push(b1);
        }
        let mut frame_detected = 0;
        for (i, carrier) in TMCC_CARRIER[segment].iter().enumerate() {
            tmcc[i].push(demap_bpsk(segment_carriers[*carrier]));
            if tmcc[i].len() == frame_size {
                let demod = demodulate_dbpsk(&tmcc[i]);
                if (demod[0] == 0b10101100 && demod[1] == 0b01110111)
                    || (demod[0] == 0b01010011 && demod[1] == 0b10001000)
                {
                    let p1 = calc_272_190_parity(&demod, 16, 120);
                    let p2 = get_bits128(&demod, 121, 202);
                    if p1 != p2 {
                        // TODO: calc syndrome
                        eprintln!("tmcc error");
                        eprintln!("{:082b}", p1);
                        eprintln!("{:082b}", p2);
                    } else {
                        frame_detected += 1;
                        let tmcc = decode_tmcc(&demod);
                        let prev_layer_a = layer_a;
                        layer_a = tmcc.current.layer_a;
                        if prev_layer_a != layer_a {
                            if prev_layer_a.carrier_modulation != layer_a.carrier_modulation {
                                panic!("{:?} is not supported", layer_a.carrier_modulation);
                            }
                            if prev_layer_a.number_of_segments != layer_a.number_of_segments {
                                panic!(
                                    "number_of_segments: {:?} is not supported",
                                    layer_a.number_of_segments
                                );
                            }
                            if prev_layer_a.time_interleaving_length
                                != layer_a.time_interleaving_length
                            {
                                eprintln!(
                                    "time_interleaving_length {:?} => {:?}",
                                    prev_layer_a.time_interleaving_length,
                                    layer_a.time_interleaving_length
                                );
                                time_deinterleaver = TimeDeinterleaver::new(
                                    mode,
                                    layer_a.time_interleaving_length.unwrap(),
                                );
                                byte_deinterleaver = ByteDeinterleaver::new();
                                bit_deinterleaver = QPSKDeinterleaver::new();
                                frame_symbols.clear();
                            }
                        }
                    }
                }
            }
        }
        if frame_detected >= TMCC_CARRIER[segment].len() / 2 {
            eprintln!("FRAME DETECTED {frame_number}");
            frame_number += 1;
            for t in &mut tmcc {
                t.clear();
            }
            if !frame_symbols.is_empty() && !frame_symbols.contains(&f32::INFINITY) {
                let depunctured = depuncture(&frame_symbols, layer_a.coding_rate.unwrap());
                let ts_chunks = depunctured.chunks(TSP_SIZE * 8);
                let mut conv_initial_state = 0;
                let mut byte_prbs_state = BYTE_PRBS_INITIAL_STATE;
                for i in ts_chunks {
                    let path_list = viterbi(i, conv_initial_state);
                    let mut buffer = Vec::<u8>::with_capacity(TSP_SIZE + 1);
                    let mut buffer_valid = true;
                    buffer.push(TS_SYNC_BYTE);
                    if path_list[TSP_SIZE - 1] != TS_SYNC_BYTE {
                        eprintln!("NOT SYNC {:08b}", path_list[TSP_SIZE - 1]);
                    }
                    if cnr_counter == 0 {
                        eprintln!("CNR = {:.02} dB", estimate_cnr(i, &convolve(&path_list)));
                    }
                    cnr_counter += 1;
                    if cnr_counter == cnr_interval {
                        cnr_counter = 0;
                    }
                    for (i, b) in path_list.iter().enumerate() {
                        byte_deinterleaver.push(i, *b);
                        if let Some(mut b) = byte_deinterleaver.pop(i) {
                            for bi in [128, 64, 32, 16, 8, 4, 2, 1] {
                                byte_prbs_state = byte_prbs(byte_prbs_state);
                                if byte_prbs_state & 0b100000000000000 != 0 && i != 203 {
                                    b ^= bi;
                                }
                            }
                            buffer.push(b);
                        } else {
                            for _ in [128, 64, 32, 16, 8, 4, 2, 1] {
                                byte_prbs_state = byte_prbs(byte_prbs_state);
                            }
                            buffer_valid = false;
                        }
                    }
                    if !buffer_valid {
                        eprintln!("byte buffer not filled");
                    } else {
                        match Decoder::new(TS_PARITY_SIZE).correct(&buffer[0..TSP_SIZE], None) {
                            Ok(buf) => {
                                if buf.data() != &buffer[0..TS_SIZE] {
                                    eprintln!("corrected");
                                }
                                writer.write_all(buf.data())?;
                            }
                            Err(e) => {
                                eprintln!("{:?} {}", e, buffer.len());
                                writer.write_all(&buffer[0..TS_SIZE])?;
                            }
                        }
                    }
                    conv_initial_state = TS_SYNC_BYTE.reverse_bits() >> 1;
                }
            } else {
                eprintln!("symbol buffer not filled");
            }
            frame_symbols.clear();
        } else {
            if tmcc[0].len() == frame_size {
                eprintln!("FRAME NOT DETECTED");
                for t in &mut tmcc {
                    t.copy_within(1.., 0);
                    t.pop();
                }
                frame_symbols.copy_within(MODE3_CARRIER_RANDOMIZE.len() * 2.., 0);
                frame_symbols.truncate(frame_symbols.len() - MODE3_CARRIER_RANDOMIZE.len() * 2);
            }
        }
        iq_buffer.copy_within(next_off.., 0);
        let next = iq_buffer.len() - next_off;
        offset += next;
        read_iq(&mut reader, &mut buffer, &mut iq_buffer[next..], args.bit)?;
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
enum SystemIdentification {
    Television,
    Sound,
    Reserved(u8),
}

#[derive(Debug, Clone, PartialEq, Eq)]
enum CarrierModulation {
    DQPSK,
    QPSK,
    QAM16,
    QAM64,
}

#[derive(Debug, Clone, PartialEq, Eq)]
enum TMCCData<T> {
    Data(T),
    Reserved(u8),
    Unused,
}

impl<T: Copy> TMCCData<T> {
    fn unwrap(&self) -> T {
        return match self {
            Self::Data(n) => *n,
            _ => panic!(),
        };
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct LayerParameter {
    carrier_modulation: TMCCData<CarrierModulation>,
    coding_rate: TMCCData<CodingRate>,
    time_interleaving_length: TMCCData<u8>,
    number_of_segments: TMCCData<u8>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct TransmissionParameter {
    partial_reception: bool,
    layer_a: LayerParameter,
    layer_b: LayerParameter,
    layer_c: LayerParameter,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct TMCC {
    system_idenfication: SystemIdentification,
    switcihg_indicator: u8,
    startup_control: bool,
    current: TransmissionParameter,
    next: TransmissionParameter,
}

fn decode_tmcc_layer(tmcc: &[u8], off: usize) -> LayerParameter {
    let carrier_modulation = match get_bits(tmcc, off + 0, off + 2) {
        0b000 => TMCCData::Data(CarrierModulation::DQPSK),
        0b001 => TMCCData::Data(CarrierModulation::QPSK),
        0b010 => TMCCData::Data(CarrierModulation::QAM16),
        0b011 => TMCCData::Data(CarrierModulation::QAM64),
        0b111 => TMCCData::Unused,
        n => TMCCData::Reserved(n as u8),
    };
    let coding_rate = match get_bits(tmcc, off + 3, off + 5) {
        0b000 => TMCCData::Data(CodingRate::Rate1_2),
        0b001 => TMCCData::Data(CodingRate::Rate2_3),
        0b010 => TMCCData::Data(CodingRate::Rate3_4),
        0b011 => TMCCData::Data(CodingRate::Rate5_6),
        0b100 => TMCCData::Data(CodingRate::Rate7_8),
        0b111 => TMCCData::Unused,
        n => TMCCData::Reserved(n as u8),
    };
    let time_interleaving_length = match get_bits(tmcc, off + 6, off + 8) {
        n @ (0b000..=0b011) => TMCCData::Data(n as u8),
        0b111 => TMCCData::Unused,
        n => TMCCData::Reserved(n as u8),
    };
    let number_of_segments = match get_bits(tmcc, off + 9, off + 12) {
        n @ 0b0001..=0b1101 => TMCCData::Data(n as u8),
        n @ (0b0000 | 0b1110) => TMCCData::Reserved(n as u8),
        _ => TMCCData::Unused,
    };
    return LayerParameter {
        carrier_modulation,
        coding_rate,
        time_interleaving_length,
        number_of_segments,
    };
}

fn decode_tmcc_parameter(tmcc: &[u8], off: usize) -> TransmissionParameter {
    let partial_reception = get_bits(tmcc, off, off) != 0;
    let layer_a = decode_tmcc_layer(tmcc, off + 1);
    let layer_b = decode_tmcc_layer(tmcc, off + 1 + 13);
    let layer_c = decode_tmcc_layer(tmcc, off + 1 + 13 + 13);
    return TransmissionParameter {
        partial_reception,
        layer_a,
        layer_b,
        layer_c,
    };
}

fn decode_tmcc(tmcc: &[u8]) -> TMCC {
    let system_idenfication = match get_bits(tmcc, 19, 20) {
        0b00 => SystemIdentification::Television,
        0b01 => SystemIdentification::Sound,
        n => SystemIdentification::Reserved(n as u8),
    };
    let switcihg_indicator = get_bits(tmcc, 21, 24) as u8;
    let startup_control = get_bits(tmcc, 25, 25) != 0;
    let current = decode_tmcc_parameter(tmcc, 26);
    let next = decode_tmcc_parameter(tmcc, 66);
    return TMCC {
        system_idenfication,
        switcihg_indicator,
        startup_control,
        current,
        next,
    };
}

fn get_bits(bytes: &[u8], begin: usize, end: usize) -> u32 {
    let mut result = 0u32;
    let mut b = 1 << (end - begin);
    for in_bi in begin..=end {
        if bytes[in_bi / 8] & (1 << (in_bi % 8)) != 0 {
            result |= b;
        }
        b >>= 1;
    }
    return result;
}

fn get_bits128(bytes: &[u8], begin: usize, end: usize) -> u128 {
    let mut result = 0u128;
    let mut b = 1 << (end - begin);
    for in_bi in begin..=end {
        if bytes[in_bi / 8] & (1 << (in_bi % 8)) != 0 {
            result |= b;
        }
        b >>= 1;
    }
    return result;
}

fn calc_272_190_parity(data: &[u8], begin: usize, end: usize) -> u128 {
    const POLY: u128 = (1 << 77)
        | (1 << 76)
        | (1 << 71)
        | (1 << 67)
        | (1 << 66)
        | (1 << 56)
        | (1 << 52)
        | (1 << 48)
        | (1 << 40)
        | (1 << 36)
        | (1 << 34)
        | (1 << 24)
        | (1 << 22)
        | (1 << 18)
        | (1 << 10)
        | (1 << 4)
        | 1;
    let mut d = 0u128;
    for in_bi in begin..=end {
        d <<= 1;
        if (data[in_bi / 8] & (1 << (in_bi % 8)) != 0) ^ ((d & (1 << 82)) != 0) {
            d ^= POLY;
        }
    }
    return d & ((1 << 82) - 1);
}
