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
package krigings;

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
import org.jgrasstools.gears.libs.exceptions.ModelsIllegalargumentException;
import org.jgrasstools.gears.libs.exceptions.ModelsRuntimeException;
import org.jgrasstools.gears.libs.modules.JGTConstants;
import org.jgrasstools.gears.libs.modules.JGTModel;
import org.jgrasstools.gears.libs.modules.ModelsEngine;
import org.jgrasstools.gears.libs.monitor.IJGTProgressMonitor;
import org.jgrasstools.gears.libs.monitor.LogProgressMonitor;
import org.jgrasstools.gears.utils.math.matrixes.ColumnVector;
import org.jgrasstools.gears.utils.math.matrixes.LinearSystem;
import org.jgrasstools.gears.utils.sorting.QuickSortAlgorithm;
import org.jgrasstools.hortonmachine.i18n.HortonMessageHandler;
import org.opengis.feature.simple.SimpleFeature;

import VGM.VGM;

import com.vividsolutions.jts.geom.Coordinate;
import com.vividsolutions.jts.geom.Geometry;

@Description("Ordinary kriging algorithm.")
@Documentation("Kriging.html")
@Author(name = "Giuseppe Formetta, Daniele Andreis, Silvia Franceschi, Andrea Antonello", contact = "http://www.hydrologis.com,  http://www.ing.unitn.it/dica/hp/?user=rigon")
@Keywords("Kriging, Hydrology")
@Label(JGTConstants.STATISTICS)
@Name("kriging")
@Status(Status.EXPERIMENTAL)
@License("General Public License Version 3 (GPLv3)")
@SuppressWarnings("nls")
public class Krigings extends JGTModel {

	@Description("The vector of the measurement point, containing the position of the stations.")
	@In
	public SimpleFeatureCollection inStations = null;

	@Description("The field of the vector of stations, defining the id.")
	@In
	public String fStationsid = null;

	@Description("The field of the vector of stations, defining the elevation.")
	@In
	public String fStationsZ = null;

	@Description("The type of theoretical semivariogram: exponential, gaussian, spherical, pentaspherical"
			+ "linear, circular, bessel, periodic, hole, logaritmic, power, spline")
	@In
	public String pSemivariogramType = null;

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


	@Description("The semivariogram mode: 0=Integral scale; 1=classical semivariogram (sill range nugget).")
	@In
	public int defaultVariogramMode = 0;

	@Description("The integral scale.")
	@In
	public double[] pIntegralscale = null;

	@Description("The variance.")
	@In
	public double pVariance = 0;

	@Description("Switch for logaritmic run selection.")
	@In
	public boolean doLogarithmic = false;

	@Description("The progress monitor.")
	@In
	public IJGTProgressMonitor pm = new LogProgressMonitor();

	@Description("Include zeros in computations (default is true).")
	@In
	public boolean doIncludezero = true;

	@Description("The range if the models runs with the gaussian variogram.")
	@In
	public double pA;

	@Description("The sill if the models runs with the gaussian variogram.")
	@In
	public double pS;

	@Description("Is the nugget if the models runs with the gaussian variogram.")
	@In
	public double pNug;

	@Description("In the case of kriging with neighbor, maxdist is the maximum distance "
			+ "within the algorithm has to consider the stations")
	@In
	public double maxdist;

	@Description("In the case of kriging with neighbor, inNumCloserStations is the number "
			+ "of stations the algorithm has to consider")
	@In
	public int inNumCloserStations;


	@Description("The hashmap withe the interpolated results")
	@Out
	public HashMap<Integer, double[]> outData = null;

	private static final double TOLL = 1.0d * 10E-8;

	private HortonMessageHandler msg = HortonMessageHandler.getInstance();




	/**
	 * Executing ordinary kriging.
	 * <p>
	 * <li>Verify if the parameters are correct.
	 * <li>Calculating the matrix of the covariance (a).
	 * <li>For each point to interpolated, evalutate the know term vector (b)
	 * and solve the system (a x)=b where x is the weight.
	 * </p>
	 * 
	 * @throws SchemaException
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
				if (defaultVariogramMode == 0) {
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
				} else if (defaultVariogramMode == 1) {
					if (doIncludezero) {
						if (Math.abs(h[0]) >= 0.0) { // TOLL
							xStationList.add(round(coordinate.x, 3));
							yStationList.add(round(coordinate.y, 3));
							zStationList.add(z);
							hStationList.add(h[0]);
							idStationList.add(id);

							n1 = n1 + 1;
						}
					} else {
						if (Math.abs(h[0]) > 0.0) { // TOLL
							xStationList.add(round(coordinate.x, 3));
							yStationList.add(round(coordinate.y, 3));
							zStationList.add(z);
							hStationList.add(h[0]);
							idStationList.add(id);

							n1 = n1 + 1;
						}
					}

				}
			}
		} finally {
			stationsIter.close();
		}

		int nStaz = xStationList.size();


		/*
		 * Check if the coordinates or the values are the same for all the measurements stations.
		 * xStationVector has the dimensions of the coordinates of the measurements points 
		 * plus 1 (the station where it is going to interpolate)
		 */


		double[] xStationVector = new double[nStaz + 1];
		double[] yStationVector = new double[nStaz + 1];
		double[] zStationVector = new double[nStaz + 1];
		double[] hStationVector = new double[nStaz + 1];
		int[] idStationVector = new int[nStaz + 1];

		boolean areAllEquals = true;
		if (nStaz != 0) {
			xStationVector[0] = xStationList.get(0);
			yStationVector[0] = yStationList.get(0);
			zStationVector[0] = zStationList.get(0);
			hStationVector[0] = hStationList.get(0);
			idStationVector[0] = idStationList.get(0);
			double previousValue = hStationVector[0];

			/*for each station added to the vector, it checks the coordinate/values and if they are
			 * not null or equal, it adds to the list of the availble stations. If the coordinates
			 * or the values are all different, the flag areAllEquals  becomes false. 
			 *  
			 */
			for (int i = 0; i < nStaz; i++) {

				double xTmp = xStationList.get(i);
				double yTmp = yStationList.get(i);
				double zTmp = zStationList.get(i);
				double hTmp = hStationList.get(i);
				int idTmp = idStationList.get(i);

				boolean doubleStation = ModelsEngine.verifyDoubleStation( xStationVector, yStationVector, zStationVector, hStationVector,
						xTmp,yTmp, zTmp, hTmp, i, false, pm);
				if (!doubleStation) {
					xStationVector[i] = xTmp;
					yStationVector[i] = yTmp;
					zStationVector[i] = zTmp;
					hStationVector[i] = hTmp;
					idStationVector[i] = idTmp;
					if (areAllEquals && hStationVector[i] != previousValue) {
						areAllEquals = false;
					}
					previousValue = hStationVector[i];
				}
			}
		}

		LinkedHashMap<Integer, Coordinate> pointsToInterpolateId2Coordinates = null;

		int numPointToInterpolate = 0;


		pointsToInterpolateId2Coordinates = getCoordinate(
				numPointToInterpolate, inInterpolate, fInterpolateid);


		Set<Integer> pointsToInterpolateIdSet = pointsToInterpolateId2Coordinates
				.keySet();
		Iterator<Integer> idIterator = pointsToInterpolateIdSet.iterator();

		int j = 0;

		int[] idArray = new int[pointsToInterpolateId2Coordinates.size()];

		double[] result = new double[pointsToInterpolateId2Coordinates.size()];

		while (idIterator.hasNext()) {
			n1 = xStationVector.length - 1;
			double sum = 0.;
			int id = idIterator.next();
			idArray[j] = id;
			Coordinate coordinate = (Coordinate) pointsToInterpolateId2Coordinates.get(id);

			// coordinate of the were it is going to interpolate
			xStationVector[n1] = coordinate.x;
			yStationVector[n1] = coordinate.y;
			zStationVector[n1] = coordinate.z;

			double[] xStation = xStationVector;
			double[] yStation = yStationVector;
			double[] zStation = zStationVector;
			double[] hStation = hStationVector;			
			int[] idStation = idStationVector;


			// in case of kriging with neighbour computes the distances between the
			// point where is going to interpolate and the other stations and it
			// sorts them
			if (inNumCloserStations > 0 || maxdist>0) {

				inNumCloserStations= (inNumCloserStations> nStaz)? nStaz:inNumCloserStations;	

				double x2, y2;
				double dDifX, dDifY;
				int iCount = xStationVector.length;
				double distanceVector[] = new double[iCount];
				double pos[] = new double[iCount];


				for (int jj = 0; jj < iCount; jj++) {

					x2 = xStationVector[jj];
					y2 = yStationVector[jj];

					dDifX = xStationVector[n1] - x2;
					dDifY = yStationVector[n1] - y2;
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
						if (Math.abs(hStationVector[(int) pos[i]]) >= 0.0) { // TOLL
							xStationWithNeighbour[n1] = xStationVector[(int) pos[i]];
							yStationWithNeighbour[n1] = yStationVector[(int) pos[i]];
							zStationWithNeighbour[n1] = zStationVector[(int) pos[i]];
							idStationWithNeighbour[n1] = idStationVector[(int) pos[i]];

							hWithNeighbour[n1] = hStationVector[(int) pos[i]];
							n1 += 1;
						}
					} else {
						if (Math.abs(hStationVector[(int) pos[i]]) > 0.0) {
							xStationWithNeighbour[n1] = xStationVector[(int) pos[i]];
							yStationWithNeighbour[n1] = yStationVector[(int) pos[i]];
							zStationWithNeighbour[n1] = zStationVector[(int) pos[i]];
							idStationWithNeighbour[n1] = idStationVector[(int) pos[i]];

							hWithNeighbour[n1] = hStationVector[(int) pos[i]];
							n1 += 1;

						}
					}

				}


				// final stations after the neighbour 
				yStationWithNeighbour[dim-1] = coordinate.y;
				xStationWithNeighbour[dim-1] = coordinate.x;
				zStationWithNeighbour[dim-1] = coordinate.z;
				idStationWithNeighbour[dim-1] = id;
				xStation = xStationWithNeighbour;
				yStation = yStationWithNeighbour;
				zStation = zStationWithNeighbour;
				idStation = idStationWithNeighbour;
				hStation = hWithNeighbour;

			}

			if (n1 != 0) {
				if (doLogarithmic) {
					for (int i = 0; i < nStaz; i++) {
						if (hStation[i] > 0.0) {
							hStation[i] = Math.log(hStation[i]);
						}
					}

				}


				if (!areAllEquals && n1 > 1) {
					pm.beginTask(msg.message("kriging.working"),
							pointsToInterpolateId2Coordinates.size());

					double h0 = 0.0;


					/*
					 * calculating the covariance matrix.
					 */
					double[][] covarianceMatrix = covMatrixCalculating(xStation, yStation, zStation, n1);

					double[] knownTerm = knownTermsCalculation(xStation,yStation, zStation, n1);

					/*
					 * solve the linear system, where the result is the weight.
					 */
					ColumnVector knownTermColumn = new ColumnVector(knownTerm);

					LinearSystem linearSystem = new LinearSystem(covarianceMatrix);

					ColumnVector solution = linearSystem.solve(knownTermColumn,true);

					double[] moltiplicativeFactor = solution.copyValues1D();

					for (int k = 0; k < n1; k++) {
						h0 = h0 + moltiplicativeFactor[k] * hStation[k];
						sum = sum + moltiplicativeFactor[k];

					}

					if (doLogarithmic) {
						h0 = Math.exp(h0);
					}
					result[j] = h0;
					j++;

					if (Math.abs(sum - 1) >= TOLL) {
						throw new ModelsRuntimeException(
								"Error in the coffeicients calculation", this
								.getClass().getSimpleName());
					}
					pm.worked(1);
				} else if (n1 == 1 || areAllEquals) {

					double tmp = hStation[0];
					pm.message(msg.message("kriging.setequalsvalue"));
					pm.beginTask(msg.message("kriging.working"),
							pointsToInterpolateId2Coordinates.size());
					result[j] = tmp;
					j++;
					n1 = 0;
					pm.worked(1);

				}

				pm.done();

			} else {

				pm.errorMessage("No value for this time step");
				j = 0;
				double[] value = inData.values().iterator().next();
				result[j] = value[0];
				j++;

			}
		}

		storeResult(result, idArray);

	}

	/**
	 * Verify the input of the model.
	 */
	private void verifyInput() {
		if (inData == null || inStations == null) {
			throw new NullPointerException( msg.message("kriging.stationProblem"));
		}

		if (defaultVariogramMode != 0 && defaultVariogramMode != 1) {
			throw new IllegalArgumentException( msg.message("kriging.variogramMode"));
		}
		if (defaultVariogramMode == 0) {
			if (pVariance == 0 || pIntegralscale[0] == 0
					|| pIntegralscale[1] == 0 || pIntegralscale[2] == 0) {

				pm.errorMessage(msg.message("kriging.noParam"));
				pm.errorMessage("varianza " + pVariance);
				pm.errorMessage("Integral scale x " + pIntegralscale[0]);
				pm.errorMessage("Integral scale y " + pIntegralscale[1]);
				pm.errorMessage("Integral scale z " + pIntegralscale[2]);
			}
		}
		if (defaultVariogramMode == 1) {
			if (pNug == 0 || pS == 0 || pA == 0) {
				pm.errorMessage(msg.message("kriging.noParam"));
				pm.errorMessage("Nugget " + pNug);
				pm.errorMessage("Sill " + pS);
				pm.errorMessage("Range " + pA);
			}
		}

		if (inInterpolate == null) {
			throw new ModelsIllegalargumentException(
					msg.message("kriging.noPoint"), this);
		}


	}


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
	 * @param nStaz
	 * @param collection
	 * @throws Exception
	 *             if a field of elevation isn't the same of the collection
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
	 * 
	 * 
	 * @param x
	 *            the x coordinates.
	 * @param y
	 *            the y coordinates.
	 * @param z
	 *            the z coordinates.
	 * @param n
	 *            the number of the stations points.
	 * @return
	 */
	private double[][] covMatrixCalculating(double[] x, double[] y, double[] z, int n) {

		double[][] covarianceMatrix = new double[n + 1][n + 1];

		if (defaultVariogramMode == 0) {
			for (int j = 0; j < n; j++) {
				for (int i = 0; i <= j; i++) {
					double rx = x[i] - x[j];
					double ry = y[i] - y[j];
					double rz = z[i] - z[j];


					covarianceMatrix[j][i] = variogram(rx, ry, rz);
					covarianceMatrix[i][j] = variogram(rx, ry, rz);

				}
			}
		} else if (defaultVariogramMode == 1) {
			for (int j = 0; j < n; j++) {
				for (int i = 0; i < n; i++) {
					double rx = x[i] - x[j];
					double ry = y[i] - y[j];
					double rz = z[i] - z[j];


					covarianceMatrix[j][i] = variogram(pNug, pA, pS, rx, ry, rz);
					covarianceMatrix[i][j] = variogram(pNug, pA, pS, rx, ry, rz);

				}
			}

		}
		for (int i = 0; i < n; i++) {
			covarianceMatrix[i][n] = 1.0;
			covarianceMatrix[n][i] = 1.0;

		}
		covarianceMatrix[n][n] = 0;
		return covarianceMatrix;

	}

	/**
	 * 
	 * @param x
	 *            the x coordinates.
	 * @param y
	 *            the y coordinates.
	 * @param z
	 *            the z coordinates.
	 * @param n
	 *            the number of the stations points.
	 * @return
	 */
	private double[] knownTermsCalculation(double[] x, double[] y, double[] z,
			int n) {

		// known terms vector 
		double[] gamma = new double[n + 1];

		if (defaultVariogramMode == 0) {
			for (int i = 0; i < n; i++) {
				double rx = x[i] - x[n];
				double ry = y[i] - y[n];
				double rz = z[i] - z[n];
				gamma[i] = variogram(rx, ry, rz);
			}
		} else if (defaultVariogramMode == 1) {
			for (int i = 0; i < n; i++) {
				double rx = x[i] - x[n];
				double ry = y[i] - y[n];
				double rz = z[i] - z[n];
				gamma[i] = variogram(pNug, pA, pS, rx, ry, rz);
			}

		}
		gamma[n] = 1.0;
		return gamma;

	}



	private double variogram(double nug, double range, double sill, double rx, double ry, double rz) {
		if (isNovalue(rz)) {
			rz = 0;
		}
		double h2 = Math.sqrt(rx * rx + rz * rz + ry * ry);
		double vgmResult;

		if(h2!=0){			
			VGM vgm=new VGM();
			vgmResult=vgm.calculateVGM(pSemivariogramType,h2, sill, range, nug);
		}else {
			vgmResult=0;
		}
		return vgmResult;
	}

	/**
	 * 
	 * @param rx
	 *            x distance.
	 * @param ry
	 *            y distance.
	 * @param rz
	 *            z distance.
	 * @return
	 */
	private double variogram(double rx, double ry, double rz) {
		if (isNovalue(rz)) {
			rz = 0;
		}
		double h2 = (rx / pIntegralscale[0]) * (rx / pIntegralscale[0])
				+ (ry / pIntegralscale[1]) * (ry / pIntegralscale[1])
				+ (rz / pIntegralscale[2]) * (rz / pIntegralscale[2]);
		if (h2 < TOLL) {
			return pVariance;
		} else {
			return pVariance * Math.exp(-Math.sqrt(h2));
		}

	}

	/**
	 * Store the result in a HashMcovarianceMatrix (if the mode is 0 or 1)
	 * 
	 * @param result2
	 *            the result of the model
	 * @param id
	 *            the associated id of the calculating points.
	 * @throws SchemaException
	 * @throws SchemaException
	 */
	private void storeResult(double[] result2, int[] id) throws SchemaException {
		outData = new HashMap<Integer, double[]>();
		for (int i = 0; i < result2.length; i++) {
			outData.put(id[i], new double[] { result2[i] });
		}
	}

}
