use std::fs::{create_dir, File};
use std::io::{BufWriter, Write};
use std::path;

pub fn write_vtk(
    file: &str,
    dx: f64,
    velocity: &[f64],
    pressure: &[f64],
    cross_section_length: &[f64],
) -> std::io::Result<()> {
    let filepath = path::Path::new(file);
    if !filepath.parent().unwrap().exists() {
        create_dir(filepath.parent().unwrap())?;
    }

    // let mut f = File::create(file).expect("Unable to open vtk file");
    let file = File::create(file).expect("Unable to open vtk file");
    let mut f = BufWriter::new(file);

    let n_points = velocity.len();

    write!(f, "# vtk DataFile Version 2.0\n\n")?;
    write!(f, "ASCII\n\n")?;
    write!(f, "DATASET UNSTRUCTURED_GRID\n\n")?;
    write!(f, "POINTS {n_points} float\n\n")?;
    for i in 0..n_points {
        writeln!(
            f,
            "{:.16e} 0.0000000000000000e+00 0.0000000000000000e+00",
            i as f64 * dx
        )?;
    }
    writeln!(f)?;

    write!(f, "POINT_DATA {n_points}\n\n")?;

    writeln!(f, "VECTORS velocity float")?;
    for v in velocity {
        writeln!(
            f,
            "{:.16e} 0.0000000000000000e+00 0.0000000000000000e+00",
            v
        )?;
    }
    writeln!(f)?;

    writeln!(f, "SCALARS pressure float")?;
    writeln!(f, "LOOKUP_TABLE default")?;
    for v in pressure {
        writeln!(f, "{:.16e}", v)?;
    }
    writeln!(f)?;

    writeln!(f, "SCALARS diameter float")?;
    writeln!(f, "LOOKUP_TABLE default")?;
    for v in cross_section_length {
        writeln!(f, "{:.16e}", v)?;
    }
    writeln!(f)?;
    Ok(())
}
