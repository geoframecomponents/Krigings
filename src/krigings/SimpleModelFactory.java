package krigings;



public class SimpleModelFactory {
	
	public static Model createModel(String type,double dist, double sill, double range, double nug){
		Model model=null;

		if ("exponential".equals(model)) {
			
			model= new Exponential(dist, sill, range, nug);
			
        }
        if ("gaussian".equals(model)) {
        	
        	model=new Gaussian (dist, sill, range, nug);
        	
        }
        if ("spherical".equals(model)) {
        	
        	model=new Spherical(dist, sill, range, nug);

        }
        if ("pentaspherical".equals(model)) {
        	
        	model=new Pentaspherical(dist, sill, range, nug);

        }
        if ("linear".equals(model)) {
        	
        	model=new Linear (dist, sill, range, nug);
        }
        if ("circular".equals(model)) {
        	
        	model=new Circolar (dist, sill, range, nug);

        }
        if ("bessel".equals(model)) {
        	
        	model=new Bessel (dist, sill, range, nug);

        }
        if ("periodic".equals(model)) {

        	model=new Periodic (dist, sill, range, nug);
        }
        if ("hole".equals(model)) {
        	
        	model=new Hole (dist, sill, range, nug);

        }

        if ("power".equals(model)) {
        	
        	model = new Power(dist, sill, range, nug);

        }
        if ("spline".equals(model)) {
        	
        	model =new Spline(dist, sill, range, nug);

        }
		
		

		return model;
		
	}

}
