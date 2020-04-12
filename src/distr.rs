use std::io::{BufReader,BufRead};
use std::fs::File;
use rand;
//lfi --> level index (0 inclusive) from index
pub fn lfi(index: usize, k: usize ) -> usize {
    if index < k {
        (index as isize - k as isize + 1).abs() as usize
    } else if index > k {
        index - k
    } else {
        0
    }
}


fn vec_add(dest_vec: &mut Vec<f64>, target_vec: &Vec<f64>) {
    let mut c = 0;
    for i in target_vec {
        dest_vec[c] += *i;
        c += 1;
    }
}

fn get_intensities(k: usize) -> Vec<Vec<f64>> {
    //reads distr.txt file and gets the required intensities
    let file = File::open("distr.txt").unwrap();
    let mut intensities: Vec<Vec<f64>> = vec![];
    for line in BufReader::new(file).lines() {
        let l = line.unwrap();
        let vec: Vec<&str> = l.split(",").collect();
        let mut intensity: Vec<f64> = vec![];
        for item in &vec {
            intensity.push(item.parse().unwrap());
        }
        intensities.push(intensity);
        if intensities.len() == 3 * k as usize {
            break;
        }
    }
    assert_eq!(intensities.len(), 3 * k as usize, "Not enough data for k layers");
    intensities
}

pub fn get_distributions(k: usize, resize: usize) -> (Vec<Vec<f64>>, Vec<Vec<f64>>) {
    let intensities = get_intensities(k);

    let mut counter = 0;
    let mut lambda: Vec<Vec<f64>> = vec![];
    let mut mju: Vec<Vec<f64>> = vec![];
    for i in &intensities {
        let level = mju.len();
        match counter {
            0 => lambda.push((*i).clone().iter().map(|x| x*resize as f64).collect()),
            1 => mju.push(i.to_vec()),
            2 => vec_add(&mut mju[level -1], i),
            _ => (),
        }

        if counter < 2 {
            counter += 1;
        } else {
            counter = 0;
            mju[level - 1] = mju[level - 1].iter().map(|x| x*resize as f64).collect::<Vec<f64>>();
        }
    }
    (lambda, mju)
}

pub struct Invariant {
    sum: f64,
    distr: f64,
}



pub fn assign_invariant(lambda: &Vec<Vec<f64>>, mju: &Vec<Vec<f64>>, k: usize) -> Vec<Vec<Invariant>> {
    let mut ratio_matrix: Vec<Vec<f64>> = Vec::new();
    let distrs_len = lambda[0].len();
    let order_num = distrs_len*2; //double the size to avoid overflow of order size
    for i in 0..k*2 {
        let index = lfi(i as usize, k);
        let mut ratio_vector: Vec<f64> = vec![];
        for j in 0..order_num {
            let denom: f64;
            if j > distrs_len - 2 {
                denom = mju[index][distrs_len - 1]
            } else {
                denom = mju[index][j + 1]
            }
            let nom: f64;
            if j < distrs_len {
                nom = lambda[index][j]
            } else {
                nom = lambda[index][distrs_len - 1]
            }
            ratio_vector.push(nom/denom);
        }
        ratio_matrix.push(ratio_vector);
    }
    let mut invariant: Vec<Vec<Invariant>> = vec![];
    for i in 0..k*2 {
        let mut distribution: Vec<f64> = vec![];
        if lfi(i as usize, k) == 0 {
            distribution.push(&ratio_matrix[i][0] * 2f64);
        } else {
            distribution.push(0f64);
        }
        let mut last_product = 1f64;
        let mut sum: f64 = distribution[0];
        for j in 1..(order_num - 1) {
            last_product *= ratio_matrix[i][j-1];
            distribution.push(last_product);
            sum += last_product;
        }
        let mut invar_vector: Vec<Invariant> = vec![];
        let mut last_sum = 0f64;
        for j in 0..(order_num - 1) {
            let mut invar: Invariant = Invariant {distr: distribution[j]/sum, sum: 0f64};
            if j < 1 {
                invar.sum = invar.distr;
                last_sum = invar.sum;
                invar_vector.push(invar);
                continue;
            }
            invar.sum = invar.distr + last_sum;
            last_sum = invar.sum;
            invar_vector.push(invar);
        }
        invariant.push(invar_vector);
    }
    invariant
}

pub fn pull_from_invariant(invariant: &Vec<Vec<Invariant>>, index: usize) -> usize {
    let r: f64 = rand::random();
    let vector_size: usize = invariant[index].len();
    for j in 0..vector_size {
        if invariant[index][j].sum >= r {
            return j;
        }
    }
    return vector_size
}

pub fn sample_from_intensity(intensity: & Vec<f64>, size: usize, stf:usize) -> isize {
    let mut num_per_stf: f64;
    let max_size = intensity.len() - 1;
    if size > max_size {
        num_per_stf = intensity[max_size];
    } else {
        num_per_stf = intensity[size];
    }
    num_per_stf *= (stf as f64) / 1000 as f64;
    let r: f64 = rand::random();
    if r < num_per_stf {
        1
    } else {
        0
    }
}

pub fn get_vec_mean(vec: &Vec<f64>) -> f64 {
    let mut sum = 0f64;
    for elem in vec {
        sum += *elem;
    }
    sum/vec.len() as f64
}

pub fn get_vec_stdv(vec: &Vec<f64>) -> f64 {
    let mean = get_vec_mean(vec);
    let mut sum = 0f64;
    for elem in vec {
        let abs = (elem - mean).abs();
        sum += abs.powi(2);
    }
    (sum/(vec.len() as f64 - 1f64)).sqrt()
}