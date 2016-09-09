package krigingsTest;

import java.util.HashMap;
import java.util.Iterator;
import java.util.Set;

import krigings.Krigings;


import org.geotools.data.simple.SimpleFeatureCollection;

import org.jgrasstools.gears.io.shapefile.OmsShapefileFeatureReader;
import org.jgrasstools.gears.io.timedependent.OmsTimeSeriesIteratorReader;
import org.jgrasstools.gears.io.timedependent.OmsTimeSeriesIteratorWriter;

import org.jgrasstools.hortonmachine.utils.HMTestCase;

//
///**
// * Test the kriging model.
// * 
// * @author daniele andreis
// * 
// */

public class TestKriging extends HMTestCase {


	/**
	 * Run the kriging models.
	 *
	 * <p>
	 * This is the case which all the station have the same value.
	 * </p>
	 * @throws Exception
	 * @throws Exception
	 */
	
	public void testKriging2() throws Exception {
		OmsShapefileFeatureReader stationsReader = new OmsShapefileFeatureReader();
		stationsReader.file = "resources/Input/rainstations.shp";
		stationsReader.readFeatureCollection();
		SimpleFeatureCollection stationsFC = stationsReader.geodata;
		//
		OmsShapefileFeatureReader interpolatedPointsReader = new OmsShapefileFeatureReader();
		interpolatedPointsReader.file = "resources/Input/basins_passirio_width0.shp";
		interpolatedPointsReader.readFeatureCollection();
		SimpleFeatureCollection interpolatedPointsFC = interpolatedPointsReader.geodata;
		//
		OmsTimeSeriesIteratorReader reader = new OmsTimeSeriesIteratorReader();
		reader.file ="resources/Input/rain_test2A.csv";
		reader.idfield = "ID";
		reader.tStart = "2000-01-01 00:00";
		reader.tTimestep = 60;
		// reader.tEnd = "2000-01-01 00:00";
		reader.fileNovalue = "-9999";
		//
		reader.initProcess();
		//
		Krigings kriging = new Krigings();
		kriging.pm = pm;
		//
		kriging.inStations = stationsFC;
		kriging.fStationsid = "ID_PUNTI_M";
		//
		kriging.inInterpolate = interpolatedPointsFC;
		kriging.fInterpolateid = "netnum";
		kriging.inNumCloserStations=5;
		//
		// it doesn't execute the model with log value.
		kriging.doLogarithmic = false;
		/*
		 * Set up the model in order to use the variogram with an explicit integral scale and
        variance.
		 */
	
	
		kriging.pVariance = 0.5;
		kriging.pIntegralscale = new double[]{10000, 10000, 100};

		//
		OmsTimeSeriesIteratorWriter writer = new OmsTimeSeriesIteratorWriter();
		writer.file = "resources/Output/kriging_interpolated.csv";
		//
		writer.tStart = reader.tStart;
		writer.tTimestep = reader.tTimestep;
		//
		while( reader.doProcess ) {
			reader.nextRecord();
			HashMap<Integer, double[]> id2ValueMap = reader.outData;
			kriging.inData = id2ValueMap;
			kriging.executeKriging();
			/*
			 * Extract the result.
			 */
	
	
			HashMap<Integer, double[]> result = kriging.outData;
			Set<Integer> pointsToInterpolateResult = result.keySet();
			Iterator<Integer> iterator = pointsToInterpolateResult.iterator();
			while( iterator.hasNext() ) {
				int id = iterator.next();
				double[] actual = result.get(id);
				assertEquals(1.0, actual[0], 0);
			}
			writer.inData = result;
			writer.writeNextLine();
		}
		//
		reader.close();
		writer.close();
	}
	// /////////////////////////////////////////////////////////////////////////////////////////
	// ///////////////////////////////FINE TEST 2
	// PASSA////////////////////////////////////////////////////
	// ///////////////////////////////////////////////////////////////////////////////////////
	//
	// /////////////////////////////////////////////////////////////////////////////////////////
	// /////////////////////////////// TEST 3
	// PASSA////////////////////////////////////////////////////
	// ///////////////////////////////////////////////////////////////////////////////////////
	// /**
	// * Run the kriging models.
	// *
	// * <p>
	// * This is the case that defaultMode=0.
	// * </p>
	// * @throws Exception
	// * @throws Exception
	// */
	public void testKriging4() throws Exception {
		OmsShapefileFeatureReader stationsReader = new OmsShapefileFeatureReader();
		stationsReader.file = "resources/Input/rainstations.shp";
		stationsReader.readFeatureCollection();
		SimpleFeatureCollection stationsFC = stationsReader.geodata;
		//
		OmsShapefileFeatureReader interpolatedPointsReader = new OmsShapefileFeatureReader();
		interpolatedPointsReader.file = "resources/Input/basins_passirio_width0.shp";
		interpolatedPointsReader.readFeatureCollection();
		SimpleFeatureCollection interpolatedPointsFC = interpolatedPointsReader.geodata;
		//
		OmsTimeSeriesIteratorReader reader = new OmsTimeSeriesIteratorReader();
		reader.file = "resources/Input/rain_test.csv";
		reader.idfield = "ID";
		reader.tStart = "2000-01-01 00:00";
		reader.tTimestep = 60;
		// reader.tEnd = "2000-01-01 00:00";
		reader.fileNovalue = "-9999";
		//
		reader.initProcess();
		//
		Krigings kriging = new Krigings();
		kriging.pm = pm;
		//
		kriging.inStations = stationsFC;
		kriging.fStationsid = "ID_PUNTI_M";
		//
		kriging.inInterpolate = interpolatedPointsFC;
		kriging.fInterpolateid = "netnum";
		//
		// it doesn't execute the model with log value.
		kriging.doLogarithmic = false;
		/*
		 * Set up the model in order to use the variogram with an explicit integral scale and
        variance.
		 */
		kriging.pVariance = 3.5;
        kriging.pIntegralscale = new double[]{10000, 10000, 100};
		//      /*
		//       * Set up the model in order to run with a FeatureCollection as point to
		//       * interpolated. In this case only 2D.
		//       */



		//
		kriging.doIncludezero = false;
		OmsTimeSeriesIteratorWriter writer = new OmsTimeSeriesIteratorWriter();
		writer.file =  "resources/Output/kriging_interpolated_2.csv";
		//
		writer.tStart = reader.tStart;
		writer.tTimestep = reader.tTimestep;
		//
		while( reader.doProcess ) {
			reader.nextRecord();
			HashMap<Integer, double[]> id2ValueMap = reader.outData;
			kriging.inData = id2ValueMap;
			kriging.executeKriging();
			/*
			 * Extract the result.
			 */
			HashMap<Integer, double[]> result = kriging.outData;
			

			//
			writer.inData = result;
			writer.writeNextLine();
		}
		//
		reader.close();
		writer.close();
	}



	/**
	 * Run the kriging models.
	 *
	 * <p>
	 * This is the case which there is only one station.
	 * </p>
	 * @throws Exception
	 * @throws Exception
	 */
	public void testKriging5() throws Exception {
		OmsShapefileFeatureReader stationsReader = new OmsShapefileFeatureReader();
		stationsReader.file = "resources/Input/rainstations.shp";
		stationsReader.readFeatureCollection();
		SimpleFeatureCollection stationsFC = stationsReader.geodata;
		//
		OmsShapefileFeatureReader interpolatedPointsReader = new OmsShapefileFeatureReader();
		interpolatedPointsReader.file = "resources/Input/basins_passirio_width0.shp";
		interpolatedPointsReader.readFeatureCollection();
		SimpleFeatureCollection interpolatedPointsFC = interpolatedPointsReader.geodata;
		//
		OmsTimeSeriesIteratorReader reader = new OmsTimeSeriesIteratorReader();
		reader.file = "resources/Input/rain_test3A.csv";
		reader.idfield = "ID";
		reader.tStart = "2000-01-01 00:00";
		reader.tTimestep = 60;
		// reader.tEnd = "2000-01-01 00:00";
		reader.fileNovalue = "-9999";
		//
		reader.initProcess();
		//
		Krigings kriging = new Krigings();
		kriging.pm = pm;
		//
		kriging.inStations = stationsFC;
		kriging.fStationsid = "ID_PUNTI_M";
		//
		kriging.inInterpolate = interpolatedPointsFC;
		kriging.fInterpolateid = "netnum";
		//
		// it doesn't execute the model with log value.
		kriging.doLogarithmic = false;
		
		 // Set up the model in order to use the variogram with an explicit integral scale and
        //variance.
		 
		kriging.pVariance = 0.5;
		kriging.pIntegralscale = new double[]{10000, 10000, 100};
		
		 // Set up the model in order to run with a FeatureCollection as point to interpolated. In this
        //case only 2D.
		 
		//
		OmsTimeSeriesIteratorWriter writer = new OmsTimeSeriesIteratorWriter();
		writer.file = "resources/Output/kriging_interpolated_3.csv";
		//
		writer.tStart = reader.tStart;
		writer.tTimestep = reader.tTimestep;
		int j = 0;
		while( reader.doProcess ) {
			reader.nextRecord();
			HashMap<Integer, double[]> id2ValueMap = reader.outData;
			kriging.inData = id2ValueMap;
			kriging.executeKriging();
			
			 // Extract the result.
			 
			HashMap<Integer, double[]> result = kriging.outData;
			
			
			
			Set<Integer> pointsToInterpolateResult = result.keySet();
			Iterator<Integer> iteratorTest = pointsToInterpolateResult.iterator();
			double expected;
			if (j == 0) {
				expected = 10.0;
			} else if (j == 1) {
				expected = 15;
			} else if (j == 2) {
				expected = 1;
			} else if (j == 3) {
				expected = 2;
			} else if (j == 4) {
				expected = 2;
			} else if (j == 5) {
				expected = 0;
			} else if (j == 6) {
				expected = 0;
			} else if (j == 7) {
				expected = 23;
			} else if (j == 8) {
				expected = 50;
			} else if (j == 9) {
				expected = 70;
			} else if (j == 10) {
				expected = 30;
			} else if (j == 11) {
				expected = 10;
			} else if (j == 12) {
				expected = 2;
			} else {
				expected = 1.0;
			}
			//
			while( iteratorTest.hasNext() ) {
				int id = iteratorTest.next();
				double[] actual = result.get(id);
				//
				assertEquals(expected, actual[0], 0);
			}
			
			
			writer.inData = result;
			writer.writeNextLine();
			j++;
		}
		//
		reader.close();
		writer.close();
	}

}
