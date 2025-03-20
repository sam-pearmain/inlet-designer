#[derive(Debug)]
pub struct Contour {
    x_coords: Vec<f64>,
    y_coords: Vec<f64>,
}

impl Contour {
    fn new() -> Self {
        Contour { x_coords: Vec::new(), y_coords: Vec::new() }
    }

    fn push_coords(&mut self, x: f64, y: f64) {
        self.x_coords.push(x);
        self.y_coords.push(y);
    }

    fn plot(&self, filename: &str) {
        todo!("not implemented plot for contours just yet")
    }
}

#[derive(Debug)]
pub struct Inlet {
    contour: Contour,
}

impl Inlet {
    pub fn export_csv(&self) {
        todo!()
    }

    pub fn plot(&self, filename: &str) {
        self.contour.plot(filename);
    }
}