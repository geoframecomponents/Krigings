/* This file is part of JGrasstools (http://www.jgrasstools.org)
 * (C) HydroloGIS - www.hydrologis.com 
 * 
 * JGrasstools is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
package trendAnalysis;

import static org.jgrasstools.gears.libs.modules.JGTConstants.isNovalue;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Set;

import oms3.annotations.Author;
import oms3.annotations.Description;
import oms3.annotations.Documentation;
import oms3.annotations.Execute;
import oms3.annotations.In;
import oms3.annotations.Keywords;
import oms3.annotations.Label;
import oms3.annotations.License;
import oms3.annotations.Name;
import oms3.annotations.Out;
import oms3.annotations.Status;

import org.geotools.data.simple.SimpleFeatureCollection;
import org.geotools.feature.FeatureIterator;
import org.geotools.feature.SchemaException;
import org.jgrasstools.gears.libs.modules.JGTModel;
import org.jgrasstools.gears.libs.modules.ModelsEngine;
import org.jgrasstools.gears.libs.monitor.IJGTProgressMonitor;
import org.jgrasstools.gears.libs.monitor.LogProgressMonitor;
import org.jgrasstools.gears.utils.sorting.QuickSortAlgorithm;
import org.jgrasstools.hortonmachine.i18n.HortonMessageHandler;
import org.opengis.feature.simple.SimpleFeature;


import com.vividsolutions.jts.geom.Coordinate;
import com.vividsolutions.jts.geom.Geometry;

import flanagan.analysis.Regression;

@Description("Ordinary kriging algorithm.")
@Documentation("Kriging.html")
@Author(name = "Marialaura Bancheri, Giuseppe Formetta, Daniele Andreis, Silvia Franceschi, Andrea Antonello")
@Keywords("Kriging, Hydrology")
@Label("")
@Name("kriging")
@Status()
@License("General Public License Version 3 (GPLv3)")
@SuppressWarnings("nls")
public class TrendAnalysis extends JGTModel {



	@Description("The vector of the measurement point, containing the position of the stations.")
	@In
	public SimpleFeatureCollection inStations = null;


	@Description("The field of the vector of stations, defining the id.")
	@In
	public String fStationsid = null;


	@Description("The field of the vector of stations, defining the elevation.")
	@In
	public String fStationsZ = null;


	@Description("The file with the measured data to be interpolated.")
	@In
	public HashMap<Integer, double[]> inData = null;


	@Description("The vector of the points in which the data have to be interpolated.")
	@In
	public SimpleFeatureCollection inInterpolate = null;


	@Description("The field of the interpolated vector points, defining the id.")
	@In
	public String fInterpolateid = null;


	@Description("The field of the interpolated vector points, defining the elevation.")
	@In
	public String fPointZ = null;


	@Description("The progress monitor.")
	@In
	public IJGTProgressMonitor pm = new LogProgressMonitor();


	@Description("Include zeros in computations (default is true).")
	@In
	public boolean doIncludezero = true;


	@Description("In the case of kriging with neighbor, maxdist is the maximum distance "
			+ "within the algorithm has to consider the stations")
	@In
	public double maxdist;


	@Description("In the case of kriging with neighbor, inNumCloserStations is the number "
			+ "of stations the algorithm has to consider")
	@In
	public int inNumCloserStations;


	@Description("Switch for detrended mode.")
	@In
	@Out
	public boolean doDetrended = true;

	
	@Description("The threshold on correlation coefficient for the trend in detrendend mode.")
	@In
	public double thresholdCorrelation;


	@Description("The hashmap withe the interpolated results")
	@Out
	public HashMap<Integer, double[]> outResiduals= null;
	
	
	@Description("The hashmap withe the interpolated results")
	@Out
	public HashMap<Integer, double[]> outDetrendedValues= null;


	private HortonMessageHandler msg = HortonMessageHandler.getInstance();



	int id;

	
	@Description("The double value of the trend")
	@Out
	double trend;


	/**
	 * Executing ordinary kriging.
	 * <p>
	 * <li>Verify if the parameters are correct.
	 * <li>Calculating the matrix of the covariance (a).
	 * <li>For each point to interpolated, evalutate the know term vector (b)
	 * and solve the system (a x)=b where x is the weight.
	 * </p>
	 *
	 * @throws Exception the exception
	 */

	@Execute
	public void executeKriging() throws Exception {

		verifyInput();

		// create the arraylist containing the station with the measurements 
		List<Double> xStationList = new ArrayList<Double>();
		List<Double> yStationList = new ArrayList<Double>();
		List<Double> zStationList = new ArrayList<Double>();
		List<Double> hStationList = new ArrayList<Double>();
		List<Integer> idStationList = new ArrayList<Integer>();

		/*
		 * counter for the number of station with measured value !=0.
		 */
		int n1 = 0;


		/*
		 * Store the station coordinates and measured data in the array.
		 * Skip data for non existing stations and also skip novalues.
		 */

		FeatureIterator<SimpleFeature> stationsIter = inStations.features();
		try {
			while (stationsIter.hasNext()) {
				SimpleFeature feature = stationsIter.next();
				int id = ((Number) feature.getAttribute(fStationsid)).intValue();

				double z = 0;
				if (fStationsZ != null) {
					try {
						z = ((Number) feature.getAttribute(fStationsZ))
								.doubleValue();
					} catch (NullPointerException e) {
						pm.errorMessage(msg.message("kriging.noStationZ"));
						throw new Exception(msg.message("kriging.noStationZ"));

					}
				}
				Coordinate coordinate = ((Geometry) feature.getDefaultGeometry()).getCentroid().getCoordinate();
				double[] h = inData.get(id);
				if (h == null || isNovalue(h[0])) {

					continue;
				}
				if (doIncludezero) {
					if (Math.abs(h[0]) >= 0.0) { // TOLL
						xStationList.add(coordinate.x);
						yStationList.add(coordinate.y);
						zStationList.add(z);
						hStationList.add(h[0]);
						idStationList.add(id);
						n1 = n1 + 1;
					}

				} else {
					if (Math.abs(h[0]) > 0.0) { // TOLL
						xStationList.add(coordinate.x);
						yStationList.add(coordinate.y);
						zStationList.add(z);
						hStationList.add(h[0]);
						idStationList.add(id);

						n1 = n1 + 1;
					}
				}
			}
		} finally {
			stationsIter.close();
		}

		int nStaz = xStationList.size();


		/*
		 * Check if the coordinates or the values are the same for all the measurements stations.
		 * xStationInitialSet has the dimensions of the coordinates of the measurements points 
		 * plus 1 (the station where it is going to interpolate)
		 */


		double[] xStationInitialSet = new double[nStaz+ 1];
		double[] yStationInitialSet = new double[nStaz+ 1];
		double[] zStationInitialSet = new double[nStaz+ 1];
		double[] hStationInitialSet = new double[nStaz+ 1];
		int[] idStationInitialSet = new int[nStaz + 1];

		boolean areAllEquals = true;
		if (nStaz != 0) {
			xStationInitialSet[0] = xStationList.get(0);
			yStationInitialSet[0] = yStationList.get(0);
			zStationInitialSet[0] = zStationList.get(0);
			hStationInitialSet[0] = hStationList.get(0);
			idStationInitialSet[0] = idStationList.get(0);
			double previousValue = hStationInitialSet[0];

			/* for each station added to the vector, it checks the coordinate/values and if they are
			 * not null or equal, it adds to the list of the availble stations. If the coordinates
			 * or the values are all different, the flag areAllEquals  becomes false.   
			 */

			for (int i = 0; i < nStaz; i++) {

				double xTmp = xStationList.get(i);
				double yTmp = yStationList.get(i);
				double zTmp = zStationList.get(i);
				double hTmp = hStationList.get(i);
				int idTmp = idStationList.get(i);

				boolean doubleStation = ModelsEngine.verifyDoubleStation( xStationInitialSet, yStationInitialSet, zStationInitialSet, 
						hStationInitialSet,xTmp,yTmp, zTmp, hTmp, i, false, pm);
				if (!doubleStation) {
					xStationInitialSet[i] = xTmp;
					yStationInitialSet[i] = yTmp;
					zStationInitialSet[i] = zTmp;
					hStationInitialSet[i] = hTmp;
					idStationInitialSet[i] = idTmp;
					if (areAllEquals && hStationInitialSet[i] != previousValue) {
						areAllEquals = false;
					}
					previousValue = hStationInitialSet[i];
				}
			}
		}

		/* in case of kriging with neighbor computes the distances between the
		 * point where is going to interpolate and the other stations and it
		 * sorts them
		 */ 

		if (inNumCloserStations > 0 || maxdist>0) {

			inNumCloserStations= (inNumCloserStations> nStaz)? nStaz:inNumCloserStations;	

			double x2, y2;
			double dDifX, dDifY;
			double distanceVector[] = new double[xStationInitialSet.length];
			double pos[] = new double[xStationInitialSet.length];


			for (int jj = 0; jj < xStationInitialSet.length; jj++) {

				x2 = xStationInitialSet[jj];
				y2 = yStationInitialSet[jj];

				dDifX = xStationInitialSet[n1] - x2;
				dDifY = yStationInitialSet[n1] - y2;
				distanceVector[jj] = Math.sqrt(dDifX * dDifX + dDifY * dDifY); 
				pos[jj] = jj;					
			}

			// sorts the distances
			QuickSortAlgorithm t = new QuickSortAlgorithm(pm);
			t.sort(distanceVector, pos);

			// posDist is the number of the stations within the distance 
			int posDist = distanceVector.length;
			for (int k = 0; k < distanceVector.length; k++) {
				if (distanceVector[k] > maxdist) {
					posDist = k;
					break;
				}
			}

			// in case there are no stations within the distance, the algorithm considers
			// at least the nearest 3
			posDist=(posDist == 1)?posDist += 3:posDist;

			/*
			 * The dimension of the new vector of the station is then defined
			 * by the actual number of the station within the distance or defined 
			 * by the users
			 */
			int dim=(inNumCloserStations > 0)?inNumCloserStations+1:posDist;


			double[] xStationWithNeighbour = new double[dim];
			double[] yStationWithNeighbour = new double[dim];
			int[] idStationWithNeighbour = new int[dim];
			double[] zStationWithNeighbour = new double[dim];
			double[] hWithNeighbour = new double[dim];


			// it is necessary to actualize the counter of the stations
			n1=0;

			for (int i = 1; i < dim; i++) {					
				if (doIncludezero) {
					if (Math.abs(hStationInitialSet[(int) pos[i]]) >= 0.0) { // TOLL
						xStationWithNeighbour[n1] = xStationInitialSet[(int) pos[i]];
						yStationWithNeighbour[n1] = yStationInitialSet[(int) pos[i]];
						zStationWithNeighbour[n1] = zStationInitialSet[(int) pos[i]];
						idStationWithNeighbour[n1] = idStationInitialSet[(int) pos[i]];

						hWithNeighbour[n1] = hStationInitialSet[(int) pos[i]];
						n1 += 1;
					}
				} else {
					if (Math.abs(hStationInitialSet[(int) pos[i]]) > 0.0) {
						xStationWithNeighbour[n1] = xStationInitialSet[(int) pos[i]];
						yStationWithNeighbour[n1] = yStationInitialSet[(int) pos[i]];
						zStationWithNeighbour[n1] = zStationInitialSet[(int) pos[i]];
						idStationWithNeighbour[n1] = idStationInitialSet[(int) pos[i]];

						hWithNeighbour[n1] = hStationInitialSet[(int) pos[i]];
						n1 += 1;

					}
				}

			}

			xStationInitialSet = xStationWithNeighbour;
			yStationInitialSet = yStationWithNeighbour;
			zStationInitialSet = zStationWithNeighbour;
			idStationInitialSet = idStationWithNeighbour;
			hStationInitialSet = hWithNeighbour;

		}


		LinkedHashMap<Integer, Coordinate> pointsToInterpolateId2Coordinates = null;

		int numPointToInterpolate = 0;


		pointsToInterpolateId2Coordinates = getCoordinate(numPointToInterpolate, inInterpolate, fInterpolateid);


		Set<Integer> pointsToInterpolateIdSet = pointsToInterpolateId2Coordinates.keySet();
		
		Iterator<Integer> idIterator = pointsToInterpolateIdSet.iterator();

		int j = 0;

		double[]  hsokres = null;
		double[] hStation=null;
		int[] idArray = new int[pointsToInterpolateId2Coordinates.size()];

		while (idIterator.hasNext()) {
			n1 = xStationInitialSet.length - 1;

			id = idIterator.next();
			idArray[j] = id;

			Coordinate coordinate = (Coordinate) pointsToInterpolateId2Coordinates.get(id);


			zStationInitialSet[n1] = coordinate.z;

			double[] zStation = zStationInitialSet;
			hStation = hStationInitialSet;			


			if (n1 != 0) {

				if (!areAllEquals && n1 > 1) {
					pm.beginTask(msg.message("analysis.working"),
							pointsToInterpolateId2Coordinates.size());


					Regression r = new Regression();
					r = new Regression(zStation, hStation);
					r.linear();
					trend = r.getBestEstimates()[0] + coordinate.z
							* r.getBestEstimates()[1];

					hsokres = r.getResiduals();

					/*If there is a trend for meteorological
					 * variables and elevation and it is statistically significant 
					 * then the residuals from this linear trend
					 * are computed for each meteorological stations.
					 */

					if (r.getXYcorrCoeff() > thresholdCorrelation && trend > 0) {

						for (int i = 0; i < hStation.length - 1; i++) {

							double t = r.getBestEstimates()[0]
									+ zStation[i] * r.getBestEstimates()[1];
							hStation[i] = hStation[i] - t;

						}
					} else {
						System.out.println("The trend is not significant");
						doDetrended=false;
					}

				}

			}

		}
		
		storeResult(hsokres,hStation, idArray);
	}



	/**
	 * Verify the input of the model.
	 */
	private void verifyInput() {
		if (inData == null || inStations == null) {
			throw new NullPointerException( msg.message("kriging.stationProblem"));
		}

	}


	/**
	 * Round.
	 *
	 * @param value is the value of the variable considered 
	 * @param places the places to consider after the comma
	 * @return the double value of the variable rounded
	 */
	public static double round(double value, int places) {
		if (places < 0)
			throw new IllegalArgumentException();

		long factor = (long) Math.pow(10, places);
		value = value * factor;
		long tmp = Math.round(value);
		return (double) tmp / factor;
	}


	/**
	 * Extract the coordinate of a FeatureCollection in a HashMap with an ID as
	 * a key.
	 *
	 * @param nStaz the number of the stations
	 * @param collection is the collection of the considered points 
	 * @param idField the field containing the id of the stations 
	 * @return the coordinate of the station
	 * @throws Exception if a field of elevation isn't the same of the collection
	 */
	private LinkedHashMap<Integer, Coordinate> getCoordinate(int nStaz,
			SimpleFeatureCollection collection, String idField)
					throws Exception {
		LinkedHashMap<Integer, Coordinate> id2CoordinatesMcovarianceMatrix = new LinkedHashMap<Integer, Coordinate>();
		FeatureIterator<SimpleFeature> iterator = collection.features();
		Coordinate coordinate = null;
		try {
			while (iterator.hasNext()) {
				SimpleFeature feature = iterator.next();
				int name = ((Number) feature.getAttribute(idField)).intValue();
				coordinate = ((Geometry) feature.getDefaultGeometry())
						.getCentroid().getCoordinate();
				double z = 0;
				if (fPointZ != null) {
					try {
						z = ((Number) feature.getAttribute(fPointZ))
								.doubleValue();
					} catch (NullPointerException e) {
						pm.errorMessage(msg.message("kriging.noPointZ"));
						throw new Exception(msg.message("kriging.noPointZ"));
					}
				}
				coordinate.z = z;
				id2CoordinatesMcovarianceMatrix.put(name, coordinate);
			}
		} finally {
			iterator.close();
		}

		return id2CoordinatesMcovarianceMatrix;
	}




	/**
	 * Store the result in a HashMap
	 *
	 * @param hsokres the vector with the residuals
	 * @param hStation the vector with the values of the stations
	 * @param id the associated id of the calculating points.
	 * @throws SchemaException the schema exception
	 */
	private void storeResult(double [] hsokres, double [] hStation, int [] id ) throws SchemaException {
		outResiduals = new HashMap<Integer, double[]>();
		outDetrendedValues = new HashMap<Integer, double[]>();
		for (int i = 0; i < hsokres.length; i++) {
			outResiduals.put(id[i], new double[] { hsokres[i] });
			outDetrendedValues.put(id[i], new double[] { hStation[i] });
		}
	}

}
