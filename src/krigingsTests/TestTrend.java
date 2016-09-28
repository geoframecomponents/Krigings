package krigingsTests;

import java.util.HashMap;


import org.geotools.data.simple.SimpleFeatureCollection;
import org.jgrasstools.gears.io.shapefile.OmsShapefileFeatureReader;
import org.jgrasstools.gears.io.timedependent.OmsTimeSeriesIteratorReader;
import org.jgrasstools.gears.io.timedependent.OmsTimeSeriesIteratorWriter;
import org.jgrasstools.hortonmachine.utils.HMTestCase;

import trendAnalysis.TrendAnalysis;


public class TestTrend extends HMTestCase {


	
	public void testTrend() throws Exception {
		OmsShapefileFeatureReader stationsReader = new OmsShapefileFeatureReader();
		stationsReader.file = "resources/Input/krigings/rainstations.shp";
		stationsReader.readFeatureCollection();
		SimpleFeatureCollection stationsFC = stationsReader.geodata;
		//

		//
		OmsTimeSeriesIteratorReader reader = new OmsTimeSeriesIteratorReader();
		reader.file ="resources/Input/krigings/rain_test.csv";
		reader.idfield = "ID";
		reader.tStart = "2000-01-01 00:00";
		reader.tTimestep = 60;
		// reader.tEnd = "2000-01-01 00:00";
		reader.fileNovalue = "-9999";
		//
		reader.initProcess();
		//
		
		TrendAnalysis trend= new TrendAnalysis();
		
		trend.fStationsid="ID_PUNTI_M";
		trend.fStationsZ="QUOTA";
		trend.thresholdCorrelation=0;
		

		trend.inStations = stationsFC;
		
		
		//
		OmsTimeSeriesIteratorWriter writer = new OmsTimeSeriesIteratorWriter();
		writer.file = "resources/Output/trend/residuals.csv";
		writer.tStart = reader.tStart;
		writer.tTimestep = reader.tTimestep;
		//
		while( reader.doProcess ) {
			reader.nextRecord();
			HashMap<Integer, double[]> id2ValueMap = reader.outData;
			trend.inData = id2ValueMap;
			
			trend.execute();
		
			
			
			/*
			 * Extract the result.
			 */
	
	
			HashMap<Integer, double[]> result = trend.outResiduals;
			
			writer.inData = result;
			writer.writeNextLine();
		}
		//
		reader.close();
		writer.close();
	}
	

}
