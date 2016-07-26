package VGM;

public class Gaussian implements Model{

	double [] dist;
	double sill;
	double range;
	double nug;


	public Gaussian (double[] dist, double sill, double range, double nug){	
		this.dist=dist;
		this.sill=sill;
		this.range=range;
		this.nug=nug;		
	}



	@Override
	public double[] result() {
		int length = dist.length;
		double[] result = new double[length];
		for (int i = 0; i < length; i++) {
			double hr;
			hr = dist[i] / (range);
			if (dist[i] != 0.0) {
				result[i] = nug + sill * (1.0 - (Math.exp(-(hr * hr))));
			}
			//System.out.println(func[i]);

		}
		return result;
	}


}
