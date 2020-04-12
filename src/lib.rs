#![crate_type = "dylib"]
extern crate rand;
extern crate libc;
use std::slice;
mod distr;

struct Queue {
    price: f64,
    size: usize
}
enum Qtype{
    Ask,
    Bid,
}
enum Age {
    Current,
    Last,
}
struct Price {
    bid: f64,
    ask: f64,
    time: u64
}

struct Orderbook {
    state: Vec<Queue>,
    last_state: Vec<Queue>,
    k: usize,
    p_ref: f64,
}

fn round64(num: f64, to: usize) -> f64 {
    (num * 10f64.powf(to as f64)).round()/10f64.powf(to as f64)
}

impl Orderbook {
    fn get_best(&self) -> (f64, f64){
        let k = self.k;
        let ask: f64;
        let bid: f64;
        if self.state[k - 1].size == 0 {
            bid = self.state[k - 2].price;
        } else {
            bid = self.state[k - 1].price;
        }
        if self.state[k].size == 0 {
            ask = self.state[k + 1].price;
        } else {
            ask = self.state[k].price;
        }
        (bid, ask)
    }

    fn print(&self) {
        println!("Limit Order-book state:");
        for q in &self.state {
            println!{"Price: {}, \t Size: {}", round64(q.price, 5), q.size};
        }
        println!("________________________");
    }
    fn get_size(&self, level: usize, age: Age, q_type: Qtype) -> usize{
        match age {
            Age::Current => match q_type {
                            Qtype::Ask => self.state[level+self.k].size,
                            Qtype::Bid => self.state[self.k-level-1].size,
                            }
            Age::Last    => match q_type {
                            Qtype::Ask => self.last_state[level+self.k].size,
                            Qtype::Bid => self.last_state[self.k-level-1].size,
                            }
            }
    }

    fn shift(&mut self, invariant: &Vec<Vec<distr::Invariant>>, tick_size: f64, direction : i32) {
        for i in 0..self.k*2 {
            self.last_state[i].price = self.state[i].price;
            self.last_state[i].size = self.state[i].size;
        }
        for i in 0..self.k*2 {
            self.state[i].price += direction as f64 *tick_size as f64;
            if i as i32 + direction >= 0 && i as i32 + direction < self.k as i32 *2i32 {
                let index: isize = i as isize + direction as isize;
                self.state[i].size = self.last_state[index as usize].size;
            }
        }
        if direction > 0 {
            self.state[self.k*2-1] = Queue {price: self.state[self.k - 2].price as f64 + tick_size,
                           size: distr::pull_from_invariant(invariant, self.k*2-1)};
        } else {
             self.state[0] = Queue {price: self.state[1].price as f64 - tick_size,
                           size: distr::pull_from_invariant(invariant, 0)}       ;
        }
        self.p_ref += direction as f64 * tick_size as f64;
        assert_eq!(self.state.len(),6);
    }

    fn get_dep_theta(&mut self, theta: f64) -> f64 {
        let mut bid_sum: f64 = 0.0;
        let mut ask_sum: f64 = 0.0;
        for i in 0..self.k {
            bid_sum += self.state[i].size as f64;
        }
        for i in self.k..self.k*2 {
            ask_sum += self.state[i].size as f64;
        }
        let ratio: f64 = ask_sum/(ask_sum + bid_sum);
        0f64.max(theta-0.0+ratio*0.0).min(1f64)
    }
}



struct Candle {
    high: f64,
    open: f64,
    close: f64,
    low: f64,
    time: u64,
    last_price_index: usize
}

struct Prices {
    prices: Vec<Price>,
    candles: Vec<Candle>
}

impl Prices {
    fn set(&mut self, lob: &Orderbook, time: u64) {
        let (bid, ask) = lob.get_best();
        (*self).prices.push(Price {bid, ask, time });
    }

    fn get_mid_price(&mut self, i: usize) -> f64 {
        (self.prices[i].ask + self.prices[i].bid)/2f64
    }

    fn add_candle(&mut self, time: u64, candle_size: usize) {
        let first_price_index: usize;
        let candle_num = self.candles.len();
        if candle_num == 0 {
            first_price_index = 0;
        } else {
            first_price_index = self.candles[candle_num - 1].last_price_index;
        }
        let last_price_index: usize = self.prices.len();
        let (open, close) = (self.get_mid_price(first_price_index),
                             self.get_mid_price(last_price_index - 1));
        let (mut high, mut low): (f64, f64) = (open, open);
        for i in first_price_index..last_price_index {
            let price = self.get_mid_price(i);
            if price > high {
                high = price;
            }
            if price < low {
                low = price;
            }
        }
        let time: u64 = (time as i64 - candle_size as i64 + 1i64) as u64;
        self.candles.push(Candle {high, open, close, low, time,
            last_price_index: last_price_index - 1});
    }

    fn print_candle(&mut self, index: usize) {
        let cn = &self.candles[index];
        let (o, h, c, l, t) = (cn.open, cn.high, cn.close, cn.low, cn.time);
        println!("open: {}, high: {}, low: {}, close: {}, time: {}", round64(o, 5), round64(h, 5),
                                                round64(l, 5), round64(c, 5), t);
    }

    fn get_mid_prices_by_time(&mut self, start: u64, end: u64, simulation_time: usize) -> Vec<f64> {
        assert!((start >= 0) & (end <= simulation_time as u64));
        let mut mid_prices = vec![];
        let mut i = 0;
        let len = self.prices.len();
        while i <= len - 1 {
            if (self.prices[i].time >= end) { break;    }
            if self.prices[i].time < start { i += 1; continue; }
            mid_prices.push(self.get_mid_price(i));
            i += 1;
        }
        mid_prices
    }

    fn get_10_min_cal_params(&mut self, stf: usize, simulation_time: usize) -> (f64, f64) {
        let min10 = 1000*60*10;
        assert!(simulation_time >= min10);
        let mut nc: usize = 0;
        let mut na: usize = 0;
        let (mrr, vol);
        let mut i = 2;
        let mut cur_time = self.prices[i].time;
        let len = self.prices.len();
        while (cur_time < min10 as u64) & (i <= len - 1) {
            cur_time = self.prices[i].time;
            if (self.get_mid_price(i-2) - self.get_mid_price(i-1) < 0f64) ==
               (self.get_mid_price(i-1) - self.get_mid_price(i) < 0f64) {
                nc += 1;
            } else {
                na += 1;
            }
            i += 1;
        }
        let mid_prices = self.get_mid_prices_by_time(0, min10 as u64, simulation_time);
        mrr = nc as f64 / (na as f64 * 2f64);
        vol = distr::get_vec_stdv(&mid_prices);
        (mrr, vol)
    }
}

fn init(lob: &mut Orderbook, p_ref: f64, invariant: &Vec<Vec<distr::Invariant>> , k: usize,
        tick_size: f64) {
    let p_lowest: f64 = p_ref - tick_size/2 as f64 - tick_size*(k as f64);
    for i in 0..k*2 {
        lob.state.push(Queue {price: p_lowest + i as f64*tick_size,
                           size: distr::pull_from_invariant(invariant, i)});
        lob.p_ref = p_ref;
    }
    if lob.last_state.len() == 0 {
        lob.k = k;
        for i in 0..k*2 {
            lob.last_state.push(Queue {price: lob.state[i].price,
                                       size: lob.state[i].size});
        }
    }
}

fn reinit(lob: &mut Orderbook, p_ref: f64, invariant: &Vec<Vec<distr::Invariant>> , k: usize,
        tick_size: f64) {
    let p_lowest: f64 = p_ref - tick_size/2 as f64 - tick_size*(k as f64);
    for i in 0..k*2 {
        lob.state[i] = (Queue {price: p_lowest + i as f64*tick_size,
                           size: distr::pull_from_invariant(invariant, i)});
        lob.p_ref = p_ref;
    }
}

fn iterate(lob: &mut Orderbook, stf: usize, k: usize, lambda: &Vec<Vec<f64>>, mju: &Vec<Vec<f64>>) {
    for i in 0..k*2 {
        lob.last_state[i].size = lob.state[i].size;
        let size: usize = lob.state[i].size;
        let s_lambda: isize = distr::sample_from_intensity(&lambda[distr::lfi(i as usize, k)],
                                                           size, stf);
        let s_mju: isize = distr::sample_from_intensity(&mju[distr::lfi(i as usize, k)], size,
                                                        stf);
        let proposed_change = s_lambda - s_mju;

        lob.state[i].size = std::cmp::max(0, (size as isize + proposed_change) as usize);

    }
}

fn maybe_modify(lob: &mut Orderbook, invariant: &Vec<Vec<distr::Invariant>>, tick_size: f64,
                theta: f64, theta_reinit: f64) {
    let ask_size = lob.get_size(0, Age::Current, Qtype::Ask);
    let bid_size = lob.get_size(0, Age::Current, Qtype::Bid);
    let z_ask: i32 = (ask_size == 0) as i32;
    let z_bid: i32 = (bid_size == 0) as i32;
    let c_ask: i32 = (0 != lob.get_size(0, Age::Last, Qtype::Ask)) as i32;
    let c_bid: i32 = (0 != lob.get_size(0, Age::Last, Qtype::Bid)) as i32;
    let direction: i32 = z_ask*c_ask - z_bid*c_bid;
    if direction != 0 {

        let mut r: f64 = rand::random();
        let lob_dep_theta:f64 = lob.get_dep_theta(theta);
        if r < lob_dep_theta {
            lob.shift(invariant, tick_size, direction);
            r = rand::random();
            if r < theta_reinit {
                let p_ref = lob.p_ref;
                let k = lob.k;
                reinit(lob, p_ref, invariant, k, tick_size);
            }
        }
    }
}

fn maybe_save_prices(lob: &Orderbook, prices: &mut Prices, i: u64) {
    let (bid, ask) = lob.get_best();
    if prices.prices[prices.prices.len() - 1].bid != bid ||
        prices.prices[prices.prices.len() - 1].ask != ask {
            prices.set(lob, i);
    }
}

fn maybe_save_candle(prices: &mut Prices, i: u64, candle_size: usize) {
    if (i + 1u64) % candle_size as u64 == 0 {
        prices.add_candle(i, candle_size)
    }
}


fn simulation(k: usize, candle_size: usize, resize: usize, p_ref: f64, tick_size: f64, stf: usize,
              simulation_time: usize, theta_reinit: f64, theta: f64, o: &mut [f64], h: &mut [f64], l: &mut [f64], c: &mut [f64], calib: &mut [f64]) {

    let (lambda, mju) = distr::get_distributions(k, resize);
    let invariant = distr::assign_invariant(&lambda, &mju, k);
    let mut lob = Orderbook {state: Vec::with_capacity(k*2), last_state: Vec::with_capacity(k*2),
                             k, p_ref
    };
    init(&mut lob, p_ref, &invariant, k, tick_size);
    let mut prices = Prices {prices: vec![], candles: vec![]};
    prices.set(&lob, 0u64);
    for i in 1..simulation_time/stf {
        iterate(&mut lob, stf, k, &lambda, &mju);
        maybe_modify(&mut lob, &invariant, tick_size, theta, theta_reinit);
        maybe_save_prices(&lob, &mut prices, i as u64);
        maybe_save_candle(&mut prices, i as u64, candle_size);
        assert_eq!(lob.state.len(), 6);
    }
    for i in 0..prices.candles.len() {
        o[i] = prices.candles[i].open;
        h[i] = prices.candles[i].high;
        l[i] = prices.candles[i].low;
        c[i] = prices.candles[i].close;
    }
    let (mrr, vol) = prices.get_10_min_cal_params(stf, simulation_time);
    calib[0] = mrr;
    calib[1] = vol;
}

fn get_slice_f64<'a>(t: * mut f64, len: libc::int32_t) -> &'a mut [f64] {
   unsafe {
        assert!(!t.is_null());
        slice::from_raw_parts_mut(t, len as usize)
    }
}


#[no_mangle]
pub extern fn sim(stf: i32, simulation_time: i32, candle_size: i32, tick_size: f64, k: i32, p_ref: f64, resize: i32, theta: f64, theta_reinit: f64,
                  o: *mut f64, h: *mut f64, l: *mut f64, c: *mut f64, len: libc::int32_t, calib: *mut f64){
    let o = get_slice_f64(o, len);
    let h = get_slice_f64(h, len);
    let l = get_slice_f64(l, len);
    let c = get_slice_f64(c, len);
    let calib = get_slice_f64(calib, 2);
    simulation(k as usize, candle_size as usize, resize as usize, p_ref, tick_size, stf as usize, simulation_time as usize, theta_reinit, theta, o, h, l, c, calib);
}



