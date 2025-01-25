#![allow(non_snake_case)]

const DSPSITENAMELEN: usize = 16;
const DSPTASKNAMELEN: usize = 16;

const WRDNC_MAX_CHANS: usize = 2;
const MAXVIQPERBIN: usize = 2;
const USERMASKLEN: usize = 512;

const RPIF_HYBRID_PULSE: u32 = 2;

// lazy_static
// static UNPACK_TABLE: [f32; 65536] = {
//     let mut idx = 0;
//     [(); 65536].map(|_| {
//         let res = unpack_iq_comp(idx);
//         idx += 1;
//         res
//     })
// };

#[derive(Debug, Clone, Copy)]
struct DspTaskID {
    iSweep: u16,
    iAuxNum: u16,
    sTaskName: [u8; 1 + DSPTASKNAMELEN],
    iScanType: u8,
}

impl Default for DspTaskID {
    #[inline]
    fn default() -> DspTaskID {
        unsafe { std::mem::zeroed() }
    }
}

#[derive(Debug, Clone, Copy)]
struct RvptsPulseInfo {
    iVersion: u8,
    iMajorMode: u8,
    iPolarization: u8,
    iPhaseModSeq: u8,
    taskID: DspTaskID,
    sSiteName: [u8; 1 + DSPSITENAMELEN],
    iAqMode: u8,
    iUnfoldMode: u8,
    iPWidthCode: u8,
    fPWidthUSec: [f32; WRDNC_MAX_CHANS],
    fBandWidthMHz: [f32; WRDNC_MAX_CHANS],
    fDBzCalib: [[f32; MAXVIQPERBIN]; WRDNC_MAX_CHANS],
    fNoiseCalib: [[f32; MAXVIQPERBIN]; WRDNC_MAX_CHANS],
    fBurstCalib: [[f32; MAXVIQPERBIN]; WRDNC_MAX_CHANS],
    iSampleSize: u16,
    iMeanAngleSync: u16,
    iFlags: u32,
    iPlaybackVersion: u16,
    fGdrOffset: f32,
    fXdrOffset: f32,
    fSyClkMHz: f32,
    fWavelengthCM: f32,
    fSaturationDBM: f32,
    fRangeMaskRes: f32,
    iRangeMask: [u16; USERMASKLEN],
    fNoiseDBm: [[f32; MAXVIQPERBIN]; WRDNC_MAX_CHANS],
    fNoiseStdvDB: [[f32; MAXVIQPERBIN]; WRDNC_MAX_CHANS],
    fNoiseRangeKM: f32,
    fNoisePRFHz: f32,
    iGparmLatchSts: [u16; 2],
    iGparmImmedSts: [u16; 6],
    iGparmDiagBits: [u16; 4],
    sVersionString: [u8; 12],
    iAntStatusMask: u32,
}

impl Default for RvptsPulseInfo {
    #[inline]
    fn default() -> RvptsPulseInfo {
        unsafe { std::mem::zeroed() }
    }
}

#[derive(Debug, Clone, Copy)]
#[allow(dead_code)]
struct RvptsPulseHdrRx {
    iDataOff: u32,
    fBurstMag: f32,
    iBurstArg: u16,
    iWrapIQ: u8,
}

impl Default for RvptsPulseHdrRx {
    #[inline]
    fn default() -> RvptsPulseHdrRx {
        unsafe { std::mem::zeroed() }
    }
}

#[derive(Debug, Clone, Copy)]
struct DspInu {
    iRoll: u16,
    iPitch: u16,
    iHead: u16,
    iRollV: u16,
    iPitchV: u16,
    iHeadV: u16,
    iLatitude: u32,
    iLongitude: u32,
    iHeight: u16,
    iVelEast: u16,
    iVelNorth: u16,
    iVelUp: u16,
}

impl Default for DspInu {
    #[inline]
    fn default() -> DspInu {
        unsafe { std::mem::zeroed() }
    }
}

#[derive(Debug, Clone, Copy)]
struct Uqibits {
    iLong: [u32; 2],
}

impl Default for Uqibits {
    #[inline]
    fn default() -> Uqibits {
        unsafe { std::mem::zeroed() }
    }
}

#[derive(Debug, Clone, Copy)]
struct RvptsPulseHdr {
    iVersion: u8,
    iFlags: u8,
    iMSecUTC: u16,
    iTimeUTC: u32,
    iBtimeAPI: u32,
    iSysTime: u32,
    iPrevPRT: u32,
    iNextPRT: u32,
    iSeqNum: u32,
    iAqMode: u8,
    iPolarBits: u8,
    iTxPhase: u16,
    iNanoUTC: u32,
    iAntStatus: u32,
    iPedAz: u16,
    iPedEl: u16,
    iAzV: u16,
    iElV: u16,
    iAz: u16,
    iEl: u16,
    iNumVecs: [u16; WRDNC_MAX_CHANS],
    iMaxVecs: [u16; WRDNC_MAX_CHANS],
    iVIQPerBin: u8,
    iTgBank: u8,
    iTgWave: u16,
    uiqPerm: Uqibits,
    uiqOnce: Uqibits,
    rx: [[RvptsPulseHdrRx; MAXVIQPERBIN]; WRDNC_MAX_CHANS],
    iUTags: u32,
    inu: DspInu,
}

impl Default for RvptsPulseHdr {
    #[inline]
    fn default() -> RvptsPulseHdr {
        unsafe { std::mem::zeroed() }
    }
}

use std::fmt::Debug;
use std::io::{Cursor, BufRead, Read, Seek, SeekFrom};
use std::error::Error;
use std::str::FromStr;

use memchr::memchr;

fn read_eol(reader: &mut Cursor<&[u8]>) -> usize {
    memchr('\n' as u8, reader.fill_buf().unwrap()).unwrap()
}

fn read_scalar<T>(prefix: &str, reader: &mut Cursor<&[u8]>) -> T
where 
    T: FromStr + Default + Debug,
    <T as FromStr>::Err: Debug
{
    if !reader.fill_buf().unwrap().starts_with(prefix.as_bytes()) {
        return Default::default();
    }

    reader.seek(SeekFrom::Current(prefix.len() as i64)).unwrap();

    let eol = read_eol(reader);
    let res = std::str::from_utf8(&reader.fill_buf().unwrap()[..eol])
        .unwrap()
        .parse::<T>()
        .unwrap_or(T::default());

    reader.seek(SeekFrom::Current(1 + eol as i64)).unwrap();

    res
}

fn read_string<const LEN: usize>(prefix: &str, reader: &mut Cursor<&[u8]>) -> [u8; LEN] {
    if !reader.fill_buf().unwrap().starts_with(prefix.as_bytes()) {
        return [0; LEN];
    }

    reader.seek(SeekFrom::Current(prefix.len() as i64)).unwrap();

    let eol = read_eol(reader);
    
    let mut res = [0; LEN];
    reader.read_exact(&mut res[..eol]).unwrap();

    reader.seek(SeekFrom::Current(1)).unwrap();

    res.try_into().unwrap()
}

fn read_array<T, const LEN: usize>(prefix: &str, reader: &mut Cursor<&[u8]>) -> [T; LEN]
where 
    T: Default + FromStr + Copy,
    <T as FromStr>::Err: Debug
{
    if !reader.fill_buf().unwrap().starts_with(prefix.as_bytes()) {
        return [T::default(); LEN];
    }

    reader.seek(SeekFrom::Current(prefix.len() as i64)).unwrap();
    let eol = read_eol(reader);

    let mut values = [T::default(); LEN];

    let buffer = std::str::from_utf8(&reader.fill_buf().unwrap()[..eol]).unwrap();

    for (i, str_val) in buffer.splitn(LEN, ' ').enumerate() {
        values[i] = str_val.parse::<T>().unwrap();
    }

    reader.seek(SeekFrom::Current(1 + eol as i64)).unwrap();

    values
}

fn read_pulse_info(reader: &mut Cursor<&[u8]>) -> Result<RvptsPulseInfo, Box<dyn Error>> {    
    let mut info = RvptsPulseInfo::default();
    
    let mut single_pulse = true;

    match reader.fill_buf().unwrap() {
        x if x.starts_with(b"rvp8PulseInfo start\n") => reader.seek(SeekFrom::Current(20)).unwrap(),
        x if x.starts_with(b"rvptsPulseInfo start\n") => reader.seek(SeekFrom::Current(21)).unwrap(),
        _ => return Err("Header does not match".into()),
    };

    info.iVersion         = read_scalar("iVersion=", reader);
    info.iMajorMode       = read_scalar("iMajorMode=", reader);
    info.iPolarization    = read_scalar("iPolarization=", reader);
    info.iPhaseModSeq     = read_scalar("iPhaseModSeq=", reader);
    info.taskID.iSweep    = read_scalar("taskID.iSweep=", reader);
    info.taskID.iAuxNum   = read_scalar("taskID.iAuxNum=", reader);
    info.taskID.iScanType = read_scalar("taskID.iScanType=", reader);
    info.taskID.sTaskName = read_string("taskID.sTaskName=", reader).try_into().unwrap();
    info.sSiteName        = read_string("sSiteName=", reader).try_into().unwrap();
    info.iAqMode          = read_scalar("iAqMode=", reader);
    info.iUnfoldMode      = read_scalar("iUnfoldMode=", reader);
    info.iPWidthCode      = read_scalar("iPWidthCode=", reader);
    info.fPWidthUSec[0]   = -1.0;
    info.fPWidthUSec[0]   = read_scalar("fPWidthUSec=", reader);
    if info.fPWidthUSec[0] > 0.0 {
        info.fBandWidthMHz[0] = read_scalar("fBandWidthMHz=", reader);
        info.fDBzCalib[0][0]  = read_scalar("fDBzCalib=", reader);
        info.fDBzCalib[0][1]  = read_scalar("fDBzCalibCx=", reader);
        info.fNoiseCalib[0]   = read_array("fNoiseCalib=", reader);
        info.fBurstCalib[0]   = read_array("fBurstCalib=", reader);
    } else {
        single_pulse = false;
        info.fPWidthUSec[0]   = read_scalar("fPWidthUSec[0]=", reader);
        info.fPWidthUSec[1]   = read_scalar("fPWidthUSec[1]=", reader);
        info.fBandWidthMHz[0] = read_scalar("fBandWidthMHz[0]=", reader);
        info.fBandWidthMHz[1] = read_scalar("fBandWidthMHz[1]=", reader);
        info.fDBzCalib[0]     = read_array("fDBzCalib[0]=", reader);
        info.fDBzCalib[1]     = read_array("fDBzCalib[1]=", reader);
        info.fNoiseCalib[0]   = read_array("fNoiseCalib[0]=", reader);
        info.fNoiseCalib[1]   = read_array("fNoiseCalib[1]=", reader);
        info.fBurstCalib[0]   = read_array("fBurstCalib[0]=", reader);
        info.fBurstCalib[1]   = read_array("fBurstCalib[1]=", reader);
    }
    info.iSampleSize      = read_scalar("iSampleSize=", reader);
    info.iMeanAngleSync   = read_scalar("iMeanAngleSync=", reader);
    info.iFlags           = read_scalar("iFlags=", reader);
    assert!(single_pulse || (info.iFlags & RPIF_HYBRID_PULSE) > 0);

    info.iPlaybackVersion = read_scalar("iPlaybackVersion=", reader);
    
    info.fGdrOffset       = read_scalar("fGdrOffset=", reader);
    info.fXdrOffset       = read_scalar("fXdrOffset=", reader);
    info.fSyClkMHz        = read_scalar("fAqClkMHz=", reader);
    info.fSyClkMHz        = read_scalar("fSyClkMHz=", reader);
    
    info.fWavelengthCM    = read_scalar("fWavelengthCM=", reader);
    info.fSaturationDBM   = read_scalar("fSaturationDBM=", reader);
    info.fRangeMaskRes    = read_scalar("fRangeMaskRes=", reader);
    info.iRangeMask       = read_array("iRangeMask=", reader);

    if single_pulse {
        info.fNoiseDBm[0]    = read_array("fNoiseDBm=", reader);
        info.fNoiseStdvDB[0] = read_array("fNoiseStdvDB=", reader);
    } else {
        info.fNoiseDBm[0]    = read_array("fNoiseDBm[0]=", reader);
        info.fNoiseDBm[1]    = read_array("fNoiseDBm[1]=", reader);
        info.fNoiseStdvDB[1] = read_array("fNoiseStdvDB[0]=", reader);
        info.fNoiseStdvDB[1] = read_array("fNoiseStdvDB[1]=", reader);
    }

    info.fNoiseRangeKM  = read_scalar("fNoiseRangeKM=", reader);
    info.fNoisePRFHz    = read_scalar("fNoisePRFHz=", reader);
    info.iGparmLatchSts = read_array("iGparmLatchSts=", reader);
    info.iGparmImmedSts = read_array("iGparmImmedSts=", reader);
    info.iGparmDiagBits = read_array("iGparmDiagBits=", reader);
    info.sVersionString = read_string("sVersionString=", reader);
    info.iAntStatusMask = read_scalar("iAntStatusMask=", reader);

    match reader.fill_buf().unwrap() {
        x if x.starts_with(b"rvp8PulseInfo end\n") => reader.seek(SeekFrom::Current(18)).unwrap(),
        x if x.starts_with(b"rvptsPulseInfo end\n") => reader.seek(SeekFrom::Current(19)).unwrap(),
        x if x.starts_with(b"rvptsPulseInfo end \n") => reader.seek(SeekFrom::Current(20)).unwrap(),
        _ => {
            println!("{:?}", reader.seek(SeekFrom::Current(0)).unwrap());
            return Err("Footer does not match".into());
        },
    };

    Ok(info)
}

fn read_pulse_hdr(reader: &mut Cursor<&[u8]>) -> Result<RvptsPulseHdr, Box<dyn Error>> {
    let mut hdr = RvptsPulseHdr::default();

    loop {
        match reader.fill_buf().unwrap() {
            x if x.starts_with(b"rvp8PulseHdr start\n") => reader.seek(SeekFrom::Current(19)).unwrap(),
            x if x.starts_with(b"rvptsPulseHdr start\n") => reader.seek(SeekFrom::Current(20)).unwrap(),
            x => {
                reader.seek_relative(1).unwrap();
                continue;
            }
        };

        break
    }

    hdr.iVersion   = read_scalar("iVersion=", reader);
    hdr.iFlags     = read_scalar("iFlags=", reader);
    hdr.iMSecUTC   = read_scalar("iMSecUTC=", reader);
    hdr.iTimeUTC   = read_scalar("iTimeUTC=", reader);
    hdr.iBtimeAPI  = read_scalar("iBtimeAPI=", reader);
    hdr.iSysTime   = read_scalar("iSysTime=", reader);
    hdr.iPrevPRT   = read_scalar("iPrevPRT=", reader);
    hdr.iNextPRT   = read_scalar("iNextPRT=", reader);
    hdr.iSeqNum    = read_scalar("iSeqNum=", reader);
    hdr.iAqMode    = read_scalar("iAqMode=", reader);
    hdr.iPolarBits = read_scalar("iPolarBits=", reader);
    hdr.iTxPhase   = read_scalar("iTxPhase=", reader);
    hdr.iNanoUTC   = read_scalar("iNanoUTC=", reader);
    hdr.iAntStatus = read_scalar("iAntStatus=", reader);
    hdr.iPedAz     = read_scalar("iPedAz=", reader);
    hdr.iPedEl     = read_scalar("iPedEl=", reader);
    hdr.iAzV       = read_scalar("iAzV=", reader);
    hdr.iElV       = read_scalar("iElV=", reader);
    hdr.iAz        = read_scalar("iAz=", reader);
    hdr.iEl        = read_scalar("iEl=", reader);

    hdr.iNumVecs[0] = 0;
    hdr.iNumVecs[0] = read_scalar("iNumVecs=", reader);

    if hdr.iNumVecs[0] > 0 {
        hdr.iMaxVecs[0] = read_scalar("iMaxVecs=", reader);
        hdr.iNumVecs[1] = 0;
        hdr.iMaxVecs[1] = 0;
    } else {
        hdr.iNumVecs[0] = read_scalar("iNumVecs[0]=", reader);
        hdr.iNumVecs[1] = read_scalar("iNumVecs[1]=", reader);
        hdr.iMaxVecs[0] = read_scalar("iMaxVecs[0]=", reader);
        hdr.iMaxVecs[1] = read_scalar("iMaxVecs[1]=", reader);
    }

    hdr.iVIQPerBin    = read_scalar("iVIQPerBin=", reader);
    hdr.iTgBank       = read_scalar("iTgBank=", reader);
    hdr.iTgWave       = read_scalar("iTgWave=", reader);
    hdr.uiqPerm.iLong = read_array("uiqPerm.iLong=", reader);
    hdr.uiqOnce.iLong = read_array("uiqOnce.iLong=", reader);

    if hdr.iMaxVecs[1] > 0 {
        hdr.rx[0][0].fBurstMag = read_scalar("RX[0][0].fBurstMag=", reader);
        hdr.rx[0][0].iBurstArg = read_scalar("RX[0][0].iBurstArg=", reader);

        hdr.rx[0][1].fBurstMag = read_scalar("RX[0][1].fBurstMag=", reader);
        hdr.rx[0][1].iBurstArg = read_scalar("RX[0][1].iBurstArg=", reader);

        hdr.rx[1][0].fBurstMag = read_scalar("RX[1][0].fBurstMag=", reader);
        hdr.rx[1][0].iBurstArg = read_scalar("RX[1][0].iBurstArg=", reader);

        hdr.rx[1][1].fBurstMag = read_scalar("RX[1][1].fBurstMag=", reader);
        hdr.rx[1][1].iBurstArg = read_scalar("RX[1][1].iBurstArg=", reader);
    } else {
        hdr.rx[0][0].fBurstMag = read_scalar("RX[0].fBurstMag=", reader);
        hdr.rx[0][0].iBurstArg = read_scalar("RX[0].iBurstArg=", reader);

        hdr.rx[0][1].fBurstMag = read_scalar("RX[1].fBurstMag=", reader);
        hdr.rx[0][1].iBurstArg = read_scalar("RX[1].iBurstArg=", reader);
    }

    hdr.iUTags         = read_scalar("iUTags=", reader);
    hdr.inu.iRoll      = read_scalar("inu.iRoll=", reader);
    hdr.inu.iPitch     = read_scalar("inu.iPitch=", reader);
    hdr.inu.iHead      = read_scalar("inu.iHead=", reader);
    hdr.inu.iRollV     = read_scalar("inu.iRollV=", reader);
    hdr.inu.iPitchV    = read_scalar("inu.iPitchV=", reader);
    hdr.inu.iHeadV     = read_scalar("inu.iHeadV=", reader);
    hdr.inu.iLatitude  = read_scalar("inu.iLatitude=", reader);
    hdr.inu.iLongitude = read_scalar("inu.iLongitude=", reader);
    hdr.inu.iHeight    = read_scalar("inu.iHeight=", reader);
    hdr.inu.iVelEast   = read_scalar("inu.iVelEast=", reader);
    hdr.inu.iVelNorth  = read_scalar("inu.iVelNorth=", reader);
    hdr.inu.iVelUp     = read_scalar("inu.iVelUp=", reader);

    match reader.fill_buf().unwrap() {
        x if x.starts_with(b"rvp8PulseHdr end\n") => reader.seek(SeekFrom::Current(17)).unwrap(),
        x if x.starts_with(b"rvptsPulseHdr end \n") => reader.seek(SeekFrom::Current(19)).unwrap(),
        x if x.starts_with(b"rvptsPulseHdr end\n") => reader.seek(SeekFrom::Current(18)).unwrap(),
        _ => {
            println!("{:?}", reader.seek(SeekFrom::Current(0)).unwrap());
            return Err("Footer does not match".into());
        },
    };

    Ok(hdr)
}

fn bin_to_deg(az: u16) -> f32 {
    az as f32 * 360.0 / 65536.0
}

fn el_bin_to_deg(el: u16) -> f32 {
    match bin_to_deg(el) {
        x if x > 270.0 => x - 360.0,
        x => x
    }
}

fn unpack_iq(reader: &mut Cursor<&[u8]>, len: usize) -> Vec<[f32; 2]> {
    lazy_static::lazy_static! {
        static ref UNPACK_TABLE: Vec<f32> = (0..=65535).map(unpack_iq_comp).collect();
    }
    
    let mut iq = Vec::with_capacity(len);
    
    let mut iqu8 = Vec::new();
    iqu8.resize(len * std::mem::size_of::<u16>(), 0);
    reader.read_exact(&mut iqu8).unwrap();

    iqu8.chunks_exact(2)
        .for_each(|ch| iq.push(u16::from_le_bytes(ch.try_into().unwrap())));

    iq.chunks_exact(2).map(|codes| {
        [
            UNPACK_TABLE[codes[0] as usize],
            UNPACK_TABLE[codes[1] as usize],
        ]
    }).collect()
}

fn unpack_iq_comp(code: u16) -> f32 {
    let f_val;

    if (code & 0xF000) > 0 {
        let mut man = (code & 0x7FF) as i32;
        let exp = ((code >> 12) & 0x00F) as i32;

        if (code & 0x0800) > 0 {
            man = (man as u32 | 0xFFFFF000) as i32;
        } else {
            man = (man as u32 | 0x00000800) as i32;
        }
    
        f_val = man as f32 * (1 << exp) as u32 as f32 / 3.3554432e7;
    } else {
        f_val = ((code as i32) << 20) as f32 / 1.759218603e13;
    }

    f_val
}

#[derive(Debug, Clone)]
struct Pulse {
    hdr: RvptsPulseHdr,
    iqh: Vec<[f32; 2]>,
    iqv: Vec<[f32; 2]>,
}

fn read_pulses(reader: &mut Cursor<&[u8]>) -> Vec<Pulse> {
    let mut pulses = Vec::new();
    
    while reader.get_ref().len() - reader.position() as usize > 450 {
        let hdr = read_pulse_hdr(reader).unwrap();
        let len = hdr.iNumVecs[0] + hdr.iNumVecs[1] + 51;
        // println!("{}", reader.position());
        // println!("{}", len);
        let iqh = unpack_iq(reader, 2 * len as usize);
        let iqv = unpack_iq(reader, 2 * len as usize);

        let pulse = Pulse {
            hdr,
            iqh,
            iqv,
        };

        pulses.push(pulse);
    }

    pulses
}

const DELTA_ELEV: f32 = 0.1;
const DELTA_RANGE: f32 = 1.0;
const N_ELEV: usize = 200;
const N_RANGE: usize = 500;

type AttenTable = [[f32; N_RANGE]; N_ELEV];

fn get_atten_table(wavelength_cm: f32) -> AttenTable {
    let mut table = [[0.0; N_RANGE]; N_ELEV];
    let wavelength_corr: f32;

    if wavelength_cm >= 5.0 {
        let fraction = (wavelength_cm - 5.0) / 5.0;
        wavelength_corr = 1.2 - fraction * 0.2;
    } else {
        let fraction = (wavelength_cm - 3.3) / 1.7;
        wavelength_corr = 1.5 - fraction * 0.3;
    }

    let mut elev: f32 = 0.0;
    for ielev in 0..N_ELEV {
        let mut range: f32 = 0.0;
        
        for irange in 0..N_RANGE {
            let atten_10cm = 
                (0.4 + 3.45 * (-elev / 1.8).exp()) *
                (1.0 - (-range / (27.8 + 154.0 * (-elev / 2.2).exp())).exp());

            let atten_db = atten_10cm * wavelength_corr;
            table[ielev][irange] = atten_db;

            range += DELTA_RANGE;
        }

        elev += DELTA_ELEV;
    }

    table
}

fn get_atten(elev: f32, rangekm: f32, table: &AttenTable) -> f32 {
    let mut ielev = ((elev / DELTA_ELEV) + 0.5) as i32;
    if ielev < 0 {
        ielev = 0;
    } else if ielev > N_ELEV as i32 - 1 {
        ielev = N_ELEV as i32 - 1;
    }

    let mut irange = ((rangekm / DELTA_RANGE) + 0.5) as i32;
    if irange > N_RANGE as i32 - 1 {
        irange = N_RANGE as i32 - 1;
    } else if irange < 0 {
        irange = 0;
    }

    table[ielev as usize][irange as usize]
}

fn calc_dbz(info: &RvptsPulseInfo, mut power: f32, range_corr: f32, atten_corr: f32) -> f64 {
    power *= 10.0f32.powf(info.fSaturationDBM / 10.0);
    let noise_power = 10.0f32.powf(info.fNoiseDBm[0][0] / 10.0);

    if power > noise_power {
        let snr_hc = (power - noise_power) / noise_power;
        (10.0 * snr_hc.log10() + info.fDBzCalib[0][0] + range_corr + atten_corr) as f64
    } else {
        -9999.0
    }
}

fn calc_vel(nyquist: f32, lag1_h: [f32; 2], lag1_v: [f32; 2]) -> f64 {
    let lag1_sum = [lag1_h[0] + lag1_v[0], lag1_h[1] + lag1_v[1]];
    let arg_vel = arg(lag1_sum);
    
    ((arg_vel / std::f32::consts::PI) * nyquist) as f64 * -1.0
}

fn arg(v: [f32; 2]) -> f32 {
    if v[0] != 0.0 || v[1] != 0.0 {
        v[1].atan2(v[0])
    } else {
        0.0
    }
}

fn calc_range(info: &RvptsPulseInfo) -> (f32, f32, usize) {
    let mut bin_count = 0;
    let mut bin_start = 0;
    let mut bin_end = 0;

    for ii in 0..512usize {
        let mask = info.iRangeMask[ii];
        if mask > 0 {
            for bit in 0..16 {
                if (1 & (mask >> bit)) > 0 {
                    let bin = bit + 16 * ii;
                    if bin_count == 0 {
                        bin_start = bin;
                    }
                    bin_end = bin;
                    bin_count += 1;
                }
            }
        }
    }

    let start_range = bin_start as f32 * info.fRangeMaskRes;
    let max_range = bin_end as f32 * info.fRangeMaskRes;
    let spacing = (max_range - start_range) / (bin_count - 1) as f32;
    (start_range, spacing, bin_count)
}

fn mean_conjugate_product<A, B>(c1: A, c2: B) -> [f32; 2] 
where A: Iterator<Item = [f32; 2]> + ExactSizeIterator,
      B: Iterator<Item = [f32; 2]> + ExactSizeIterator
{
    let mut sum_i = 0.0;
    let mut sum_q = 0.0;

    let len = std::cmp::min(c1.len(), c2.len());

    for (iq1, iq2) in c1.into_iter().zip(c2.into_iter()) {
        sum_i += (iq1[0] * iq2[0]) + (iq1[1] * iq2[1]);
        sum_q += (iq1[1] * iq2[0]) - (iq1[0] * iq2[1]);
    }

    [sum_i / len as f32, sum_q / len as f32]
}

fn find_ppi_or_rhi(pulses: &Vec<Pulse>) -> bool { // true = ppi, false = rhi
    let mut ppi_diff = 0.0;
    let mut rhi_diff = 0.0;

    let mut last_az = bin_to_deg(pulses[0].hdr.iAz);
    let mut last_el = el_bin_to_deg(pulses[0].hdr.iEl);

    for pulse in pulses {
        let az = bin_to_deg(pulse.hdr.iAz);
        let el = el_bin_to_deg(pulse.hdr.iEl);

        ppi_diff += angle_diff(last_az, az).abs();
        rhi_diff += angle_diff(last_el, el).abs();

        last_az = az;
        last_el = el;
    }

    ppi_diff >= rhi_diff
}

fn angle_diff(a1: f32, a2: f32) -> f32 {
    (a2 - a1 + 180.0).rem_euclid(360.0) - 180.0
}

fn read_sweep(pulses: Vec<Pulse>, info: &RvptsPulseInfo, is_ppi: bool) -> silv::Sweep {
    // let azs = pulses.iter().enumerate().fold(BTreeMap::new(), |mut acc, (i, pulse)| {
    //     let az = U16F16::from_num(bin_to_deg(pulse.hdr.iAz).rem_euclid(360.0));
    //     acc.entry(az).or_insert(vec![i]).push(i);
    //     acc
    // });

    let get_az = |pulse: &Pulse| {
        if is_ppi {
            bin_to_deg(pulse.hdr.iAz)
        } else {
            el_bin_to_deg(pulse.hdr.iEl)
        }
    };

    let mut last_az = get_az(&pulses[0]);
    let mut dir = 0.0;
    let mut splits = vec![0];
    let mut acc_az = 0.0;
    for i in 0..pulses.len() {
        let az = get_az(&pulses[i]);
        
        let diff = angle_diff(last_az, az);
        
        // println!("{az}");
        // if  {
        // }

        let next_angle = {
            let mut angle = 0.0;

            for j in i + 1..pulses.len() {
                let a = get_az(&pulses[j]);
                if a != az {
                    angle = a;
                    break;
                }
            }
            
            angle
        };

        let diff_next = angle_diff(az, next_angle);

        if dir == 0.0 && diff != 0.0 {
            dir = diff.signum();
        } else if (dir != diff.signum() && diff != 0.0 && diff_next.signum() == diff.signum()) || diff.abs() > 10.0 || acc_az >= 360.0 {
            splits.push(i);
            dir = diff.signum();
            acc_az = 0.0;
        }

        last_az = az;
        acc_az += diff;
    }

    splits.push(pulses.len());

    for (i, sweep) in splits.windows(2).enumerate() {
        println!("Found sweep {} with {} rays starting at az: {} el: {} and ending at az: {} el: {}", 
            i, sweep[1] - sweep[0],
            bin_to_deg(pulses[sweep[0]].hdr.iAz), el_bin_to_deg(pulses[sweep[0]].hdr.iEl),
            bin_to_deg(pulses[sweep[1] - 1].hdr.iAz), el_bin_to_deg(pulses[sweep[1] - 1].hdr.iEl));
    }

    // for split in splits.iter() {
    //     for i in split.saturating_sub(10)..std::cmp::min(split + 10, pulses.len()) {
    //         println!("{}", get_az(&pulses[i]));
    //     }
    //     println!("----");
    // }

    // println!("{splits:?}");

    let (start, end) = (splits[SWEEP_N], splits[SWEEP_N + 1]);

    let (start_range, spacing, bins) = calc_range(info);
    let corr = calc_range_corr(start_range, spacing, bins);

    let elev = pulses[start..end]
        .iter()
        .map(|pulse| if is_ppi { el_bin_to_deg(pulse.hdr.iEl) } else { 0.0 })
        .sum::<f32>() / (end - start) as f32;

    let prt = pulses[start..end]
        .iter()
        .map(|pulse| pulse.hdr.iPrevPRT as f32)
        .sum::<f32>() / (end - start) as f32 / (info.fSyClkMHz * 1e6);

    let nyquist = info.fWavelengthCM / 100.0 / prt / 4.0;

    let atten_table = get_atten_table(info.fWavelengthCM);
    let atten_corr: Vec<f32> = (0..bins).map(|ii| {
        let rangekm = start_range / 1000.0 + ii as f32 * spacing / 1000.0;
        get_atten(elev, rangekm, &atten_table)
    }).collect();

    let mut sweep = silv::Sweep {
        rays: Vec::new(),
        elevation: elev,
        nyquist_velocity: nyquist,
        latitude: 35.2435,
        longitude: -97.47,
        ..Default::default()
    };

    // let azs: Vec<(U16F16, Vec<usize>)> = pulses[start..end].iter().enumerate().fold(BTreeMap::new(), |mut acc, (i, pulse)| {
    //     let az = U16F16::from_num(bin_to_deg(pulse.hdr.iAz).rem_euclid(360.0));
    //     acc.entry(az).or_insert(vec![i]).push(i);
    //     acc
    // }).into_iter().collect();

    // for chunk in azs.chunks(SAMPLES) {
    //     let num = chunk.iter().map(|(_, is)| is.len()).sum::<usize>() as f32;
    //     println!("{:?}", chunk.iter().map(|(az, _)| az).collect::<Vec<_>>());
    //     let az = chunk[0].0;
    //     // println!("{}", az);
    //     let powers = chunk.into_iter()
    //         .flat_map(|(_, is)| is)
    //         .fold(vec![0.0f32; pulses[start].iqh.len()], |mut acc, &i| {
    //             for (sum, iq) in acc.iter_mut().zip(pulses[i].iqh.iter()) {
    //                 *sum += (iq[0] * iq[0] + iq[1] * iq[1]) / num;
    //             }

    //             acc
    //     });

    //     let data = powers.iter()
    //         .zip(corr.iter().zip(atten_corr.iter()))
    //         .map(|(&power, (&corr, &atten_corr))| calc_dbz(&info, power, corr, atten_corr))
    //         .collect();

    //     let ray = silv::Ray {
    //         azimuth: az.checked_to_num().unwrap(),
    //         data: HashMap::from([("REF".into(), data)]),
    //         ..Default::default()
    //     };

    //     sweep.rays.push(ray);
    // }

    // --- //

    for chunk in pulses[start..end].chunks(SAMPLES) {
        // let az = chunk.iter().map(|pulse| bin_to_deg(pulse.hdr.iAz)).sum::<f32>() / SAMPLES_FL;
        let powers = chunk.iter()
            .fold(vec![0.0f32; pulses[start].iqh.len()], |mut acc, pulse| {
                for (sum, iq) in acc.iter_mut().zip(pulse.iqh.iter()) {
                    *sum += (iq[0] * iq[0] + iq[1] * iq[1]) / SAMPLES_FL;
                }

                acc
        });

        let dbz = powers.iter()
            .zip(corr.iter().zip(atten_corr.iter()))
            .map(|(&power, (&corr, &atten_corr))| calc_dbz(&info, power, corr, atten_corr))
            .inspect(|&v| if !v.is_finite() || v.abs() > 1e5 { println!("Found bad value") })
            .collect();

        let vel = (0..pulses[start].iqh.len()).map(|i| {
            let lag1_h = mean_conjugate_product(
                chunk.iter().skip(1).map(|c| c.iqh[i]),
                chunk.iter().map(|c| c.iqh[i])
            );

            let lag1_v = mean_conjugate_product(
                chunk.iter().skip(1).map(|c| c.iqv[i]),
                chunk.iter().map(|c| c.iqv[i])
            );

            calc_vel(nyquist, lag1_h, lag1_v)
        }).collect();

        let az = get_az(&chunk[0]);

        // if (az - 9.2).abs() < 1.0 {
        //     println!("{}", az);
        // }

        let ray = silv::Ray {
            azimuth: if is_ppi { az } else { -az + 90.0 },
            data: HashMap::from([("REF".into(), dbz), ("VEL".into(), vel)]),
            ..Default::default()
        };

        sweep.rays.push(ray);
    }

    sweep
}

fn calc_range_corr(start: f32, spacing: f32, bins: usize) -> Vec<f32> {
    let mut corr = Vec::with_capacity(bins);
    
    for ii in 0..bins {
        let range = start / 1000.0 + ii as f32 * spacing / 1000.0;
        if range < 0.001 {
            corr.push(0.0);
        } else {
            corr.push(20.0 * range.log10());
        }
    }

    corr
}

use std::collections::HashMap;

const SAMPLES: usize = 5;
const SAMPLES_FL: f32 = SAMPLES as f32;
const SWEEP_N: usize = 0; // Which sweep to output, if there are multiple

fn main() {
    let file = std::fs::read(r"C:\Users\super\Downloads\FOP1_RVP.20110524.222721.936.vcp12.1.H.460").unwrap();
    let mut reader = Cursor::new(file.as_slice());

    let info = read_pulse_info(&mut reader).unwrap();
    let pulses = read_pulses(&mut reader);

    let mut radar = silv::RadarFile {
        name: "NOP4".into(),
        sweeps: Vec::new(),
        params: HashMap::new(),
    };

    let (start, spacing, _) = calc_range(&info);

    radar.params.insert("REF".into(), silv::ParamDescription {
        description: "".into(),
        units: "".into(),
        meters_to_first_cell: start,
        meters_between_cells: spacing,
    });

    radar.params.insert("VEL".into(), silv::ParamDescription {
        description: "".into(),
        units: "".into(),
        meters_to_first_cell: start,
        meters_between_cells: spacing,
    });

    let is_ppi = find_ppi_or_rhi(&pulses);

    radar.sweeps.push(read_sweep(pulses, &info, is_ppi));

    let options = silv::RadyOptions { sort_rays_by_azimuth: true, ..Default::default() };

    silv::write(radar, ".", &options);
}

