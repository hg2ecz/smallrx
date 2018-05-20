use std::env;
use std::io::{self, Read, Write, BufReader, BufWriter};

const SAMP_RATE: i32 = 240000;
const OUTPUT_RATE: i32 = 48000;

const DECIMATE_TRANSITION_BW: f32 = 800.;
const SSB_BW: f32 = 3000.;
const AMFM_BW: f32 = 12000.;

const M_PI: f32 = std::f32::consts::PI;

#[derive(Copy,Clone)]
pub struct Cplx {
    pub re: f32,
    pub im: f32,
}

pub fn cadd(a: Cplx, b: Cplx) -> Cplx {
    let re = a.re + b.re;
    let im = a.im + b.im;
    Cplx {re, im}
}

pub fn csub(a: Cplx, b: Cplx) -> Cplx {
    let re = a.re - b.re;
    let im = a.im - b.im;
    Cplx {re, im}
}

pub fn cmul(a: Cplx, b: Cplx) -> Cplx {
    let re = a.re * b.re - a.im * b.im;
    let im = a.re * b.im + a.im * b.re;
    Cplx {re, im}
}

pub fn cabs(a: Cplx) -> f32 {
    (a.re * a.re + a.im * a.im).sqrt()
}


fn hamming(x: f32) -> f32 {
    0.54-0.46*(2.*M_PI*(0.5+(x/2.)))
}

fn main() { //syntax: rx <center_freq> <rx_freq> <l[sb]|u[sb]|a[m]|f[m]>
    let args: Vec<String> = env::args().collect();
    if args.len() != 4 {
        eprintln!("Usage: {} <center_freq> <rx_freq> <l[sb]|u[sb]|a[m]|f[m]>", args[0]);
        return;
    }
    let center_freq: f32 = args[1].parse().unwrap();
    let rx_freq: f32 = args[2].parse().unwrap();
    let modulation = args[3].chars().next().unwrap(); // l[sb], u[sb], a[m], f[m]

    let decimate_factor: usize = (SAMP_RATE / OUTPUT_RATE) as usize;
    let decimate_taps_length: usize = (4./(DECIMATE_TRANSITION_BW as f32/SAMP_RATE as f32)) as usize | 1; //(should be odd)
    let decimate_taps_middle: usize = decimate_taps_length/2;
    let dshift: f32 = 2.*std::f64::consts::PI as f32*(rx_freq - center_freq) as f32/SAMP_RATE as f32;
    eprintln!("decimate_taps_len: {}", decimate_taps_length);
    //calculate filter taps
    let mut decimate_taps: Vec<Cplx> = vec![Cplx{re: 0., im: 0.}; decimate_taps_length];
    let decimate_cutoff_rate: f32 = if modulation == 'u' || modulation == 'l' {
        (SSB_BW/2.)/SAMP_RATE as f32
    } else {
        (AMFM_BW/2.)/SAMP_RATE as f32
    };
    decimate_taps[decimate_taps_middle] = Cplx{re: 2. * M_PI * decimate_cutoff_rate * hamming(0.), im: 0.};
    for i in 1 .. decimate_taps_middle+1 {
        decimate_taps[decimate_taps_middle+i] =
           Cplx{re: (2.*M_PI*decimate_cutoff_rate*i as f32).sin()/i as f32 * hamming(i as f32/decimate_taps_middle as f32), im: 0.};
        decimate_taps[decimate_taps_middle-i] = decimate_taps[decimate_taps_middle+i]
    }
    let mut shift: f32 = 0.;
    // make the filter asymmetric in case of SSB
    let decimate_dshift: f32 = if modulation == 'u' { 1. } else { -1. } * ((SSB_BW/2.)/SAMP_RATE as f32)*2.*M_PI; // USB or LSB

    if true && (modulation == 'u' || modulation == 'l') {
        for i in 0 .. decimate_taps_length {
            decimate_taps[i] = cmul(decimate_taps[i as usize], Cplx{re: shift.sin(), im: shift.cos()});
            shift += decimate_dshift;
            if shift > 2.*M_PI { shift -= 2.*M_PI; }
        }
    }
    // normalize filter
    let mut decimate_taps_sum: f32 = 0.;
    for i in 0 .. decimate_taps_length { decimate_taps_sum += cabs(decimate_taps[i]); }
    for i in 0 .. decimate_taps_length {
        decimate_taps[i].re /= decimate_taps_sum;
        decimate_taps[i].im /= decimate_taps_sum;
    }
/*
#ifdef PRINTFREQZ
    FILE* filterfile = fopen("test.m", "w"); // <optional> <-- print the filter taps to file
    fprintf(filterfile, "#!/usr/bin/octave\n\nfreqz([\n");
    for(int i=0; i<decimate_taps_length; i++) fprintf(filterfile, " %f+(%f*i)\n", creal(decimate_taps[i]), cimag(decimate_taps[i]));
    fprintf(filterfile, "]);\ninput(\"\");\n");
    fchmod(fileno(filterfile), 0755);
    fclose(filterfile);                      // </optional>
#endif
*/
    shift = 0.;
    let mut samplebuf: Vec<Cplx> = vec![Cplx{re: 0., im: 0.}; decimate_taps_length];
    let mut last_phi: f32 = 0.;
    //let getchar: Option<i32> = std::io::stdin().bytes().next().and_then(|result| result.ok()).map(|byte| byte as i32);

    let mut input = BufReader::new(io::stdin()).bytes();
    let mut outstream = BufWriter::new(io::stdout());
    loop {
        // read 8bit[2] input from stdin   and  convert to complex float
        let in_a: u8;
        let in_b: u8;
        if let Some(in_aa) = input.next() { in_a = in_aa.unwrap() } else { break; }
        if let Some(in_bb) = input.next() { in_b = in_bb.unwrap() } else { break; }

        let sample: Cplx = Cplx {
            re: in_a as f32/(std::u8::MAX as f32/2.)-1.,
            im: in_b as f32/(std::u8::MAX as f32/2.)-1.
        };
        // oscillator & mixer
        shift += dshift;
        if shift > 2.*M_PI { shift-=2.*M_PI; }
        samplebuf.push( cmul(Cplx{re: shift.sin(), im: shift.cos()}, sample) ); // mix & buf

        // decimator: every N. sample   and   demodulate
        if samplebuf.len() >= decimate_taps_length {
            // decimate
            let mut decim: Cplx = Cplx{re: 0., im: 0.};
            for i in 0 .. decimate_taps_length { decim = cadd(decim, cmul(samplebuf[i], decimate_taps[i])); }
            samplebuf = samplebuf[decimate_factor .. decimate_taps_length].into();

            // demodulate
            let soundout: i16;
            match modulation {
                'f' => { // <-- fmdemod
                    let phi: f32 = if decim.re != 0. { (decim.im / decim.re).atan() } else { M_PI/2. };
                    let mut dphi: f32 = phi-last_phi;
                    last_phi = phi;
                    while dphi < -M_PI { dphi += 2.*M_PI; }
                    while dphi >  M_PI { dphi -= 2.*M_PI; }
                    soundout = ((std::i16::MAX-1) as f32 * (dphi/M_PI)) as i16;
                },
                'a' => soundout = (cabs(decim) * std::i16::MAX as f32) as i16, // <-- amdemod
                _   => soundout = (decim.re * std::i16::MAX as f32) as i16, // <-- ssbdemod
            }
            // write demodulated sound to stdout
            let out: [u8; 2] = [(soundout & 0xff) as u8, (soundout >> 8) as u8];
            outstream.write(&out).unwrap();
        }
    } // while(1)
}
