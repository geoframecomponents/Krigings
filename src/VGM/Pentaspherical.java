package VGM;

public class Pentaspherical implements Model{
	
	double [] dist;
	double sill;
	double range;
	double nug;
	
	
	public Pentaspherical (double[] dist, double sill, double range, double nug){	
		this.dist=dist;
		this.sill=sill;
		this.range=range;
		this.nug=nug;		
	}
	
	

	@Override
	public double[] result() {
        int length = dist.length;
        double[] result = new double[length];
        double hr = 0, h2r2;
        for (int i = 0; i < length; i++) {
            hr = dist[i] / (range);
            h2r2 = hr * hr;
            if (dist[i] != 0.0) {
                result[i] = nug + sill * (hr * ((15.0 / 8.0) + h2r2 * ((-5.0 / 4.0) + h2r2 * (3.0 / 8.0))));
            }
            if (dist[i] >= range) {
                result[i] = sill;
            }
            //System.out.println(func[i]);
        }
        return result;
	}

}
