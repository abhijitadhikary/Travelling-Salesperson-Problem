/*
 * NAME: ABHIJIT ADHIKARY
 * STUDENT ID: M00641486
 * TRAVELLING SALESPERSON PROBLEM SOLVER
 * COURSEWORK 1 FOR CSD3939 ARTIFICIAL INTELLENGENCE
 * MIDDLESEX UNIVERSITY
 * 
 * ALGORITHM REFERENCE:
 * 	NEAREST NEIGHBOR: https://en.wikipedia.org/wiki/Nearest_neighbour_algorithm
 * 	TWO OPT: https://en.wikipedia.org/wiki/2-opt
 * 	DYNAMIC PROGRAMMING: https://github.com/williamfiset/Algorithms/blob/master/com/williamfiset/algorithms/graphtheory/TspDynamicProgrammingIterative.java
 *
 * RUNNING INSTRUCTIONS:
 * 
 * 	To run a certain txt file rename the "FILE_NAME" with appropreate filename and location
 * 	
 * 	Run the desiered algorithm from the init() function and comment out the other two for better performance
 * 	Do not run nearestNeighbor() and twoOptOptimizer() at the same time from init() as nearestNeighbor() runs within twoOptOptimizer()
 * 	
 */

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

public class TSP {
	// final variables
	// test files: test1tsp.txt, test2atsp.txt, test3atsp.txt
	// final test files: TSP1-18.txt, TSP2-18.txt, TSP3-18.txt, TSP4-18.txt
	private static final String FILE_NAME = "TSP4-18.txt";
    private static final int FIRST_CITY = 0;
    private static final double MINIMUM_DISTANCE = 0;
    private static final double MAXIMUM_DISTANCE = Double.POSITIVE_INFINITY;
    private static final int NOT_VISITED = -1;
    private static final int INITIAL_IMPROVEMENT = 0;
    private static final int CITY_COLUMN = 0;
    private static final int X_COORDINATE_COLUMN = 1;
    private static final int Y_COORDINATE_COLUMN = 2;
    private static final int DATA_FIELDS = 3;
    private static final int IMPROVEMENT_LIMIT = 20;

    // manipulatable variables
    private static ArrayList<double[]> inputList;
    private static double[][] dataArray;
    private static int numCities;
    private static double[][] distanceMatrix;
    private static boolean[] visited;
    private static ArrayList<Integer> optimalRoute;
    private static long startTime;
    private static long endTime;
    private static long totalTime;
    
    private static double optimalDistance = Double.POSITIVE_INFINITY;

    public static void main(String[] args) {
    	// start the clock timer
    	startTime = System.nanoTime();
        
    	// function call to read in the data file
    	readIn();
        
    	// process the input data and runs the search algorithms
        init();

        // prints the execution time of the program
        printExecutionTime();
        
        // closes the program
        System.exit(0);
    }
    
    public static void init() {
    	
    	// function call to process the input data
    	processInputData();
        
        // function call to print the data
    	printData();
    	
    	// builds the distance matrix (city to city distance, with main diagonal as 0s)
    	buildDistanceMatrix();
    	
        // prints the distance matrix
        printDistanceMatrix();
        
        // set all the elements of "visited" array to false
        setAllToUnvisited();
        
        // runs te nearest neighbor algorithm
        //nearestNeighbor();
        
        // runs the two opt optimizer algorithm on top of the nearest neighbor algorithm
        twoOptOptimizer();
        
        // runs the dynamic programming algorithm
        dynamicProtrammingAlgorithm();

    }
    
    public static void readIn() {
        
        BufferedReader reader = null;
        
        try {
        	reader = new BufferedReader(new FileReader(FILE_NAME));
            String inputLine;
            inputList = new ArrayList<double[]>();

            
            // takes input from file until reaches the very city
            while ((inputLine = reader.readLine()) != null) {
            	// triggers if the line is not blank
                if (!inputLine.equals("")) {
                	// replaces any length of blank space with just one blank space
                	String tempString = inputLine.replaceAll("\\s{2,}", " ").trim();
                	// replaces all the tabs(\t) present in the dataset and replaces them with a white space and then trims any leading or following blank spaces
                	tempString = tempString.replaceAll("\\t", " ").trim();
                	// splits the String where blank space is present and stores in the parts String
                    String[] parts = tempString.split(" ");
                    
                    double[] tempArray = new double[parts.length];
                    
                    // converts each element of the parts String to double and stores it in the tempArray
                    for(int element = 0; element < parts.length; element++) {
                       //System.out.println(parts[element]);
                    	tempArray[element] = Double.parseDouble(parts[element]);
                    }
                    
                    // adds the tempArray to the inputList ArrayList
                    inputList.add(tempArray);
                }
            }
            
            // if the file has no data in it the program exits
            if (inputList.isEmpty()) {
            	System.out.println("\"" + FILE_NAME + "\" is Empty\nPlease enter a valid file");
            	System.exit(0);
            }

        } catch (IOException exception) {
        	// triggers if the specified file is not found
            System.out.println(exception);
            System.out.println("Could not find a file with the name \"" + FILE_NAME + "\"");
            System.out.println("Please make sure the file name is entered correctly");
            System.exit(0);

        } finally {
            try {
                if (reader != null) {
                	// closes the data file
                	reader.close();
                }
            } catch (IOException exception) {
            	// triggers if there is a problem closing the file
                System.out.println(exception);
                System.exit(0);
            }
        }
    }
    
    public static void processInputData() {
    	// assigns the numCities variable with the number of cities
        numCities = inputList.size();
        
        // initializes the visited array and the optimalRoute and the dataArray
        visited = new boolean[numCities];
        optimalRoute = new ArrayList<Integer>();
        
        dataArray = new double[numCities][DATA_FIELDS];
        
        // converts the inputList to an array
        // for each row: index 0 -> city ; 1 -> x-coordinate; 2 -> y-coordinate 
        for (int city = 0 ; city < numCities; city++) {
            dataArray[city] = inputList.get(city);
        }

    }
    
    public static void printData() {
    	// prints the input data
        System.out.println("********************************");
        System.out.println("<<<<<<<<<< INPUT DATA >>>>>>>>>>");
        System.out.println("********************************\n");
        System.out.println("CITY\tX\tY");
        System.out.println("************************");
        
        // prints for each row: index 0 -> #city ; 1 -> x-coordinate; 2 -> y-coordinate 
        // converts from double to int
        for (double[] data : dataArray) {
            System.out.print((int)data[CITY_COLUMN]);
            System.out.print("\t" + (int)data[X_COORDINATE_COLUMN]);
            System.out.println("\t" + (int)data[Y_COORDINATE_COLUMN]);
        }
    }
    
    public static double getEuclidianDistance(double x1, double y1, double x2, double y2) {
        // returns the euclidian distance between two cities
    	return Math.sqrt(Math.abs(x1-x2) * Math.abs(x1-x2) + Math.abs(y1-y2) * Math.abs(y1-y2));
    }
    
    public static void buildDistanceMatrix() {
        // builds the city to city distance matrix
        distanceMatrix = new double[numCities][numCities];
        
        for (int row = 0; row < numCities; row++) {
            for (int column = 0; column < numCities; column++) {
                if (row == column) {
                    distanceMatrix[row][column] = 0;
                } else {
                    distanceMatrix[row][column] = getEuclidianDistance(dataArray[row][X_COORDINATE_COLUMN], dataArray[row][Y_COORDINATE_COLUMN], dataArray[column][X_COORDINATE_COLUMN], dataArray[column][Y_COORDINATE_COLUMN]);
                }
            }
        }
    }
    
    public static void printDistanceMatrix() {
    	// prints the city to city distance matrix; all the elements of the main diagonal are 0s
        System.out.println("\n\n*************************************");
        System.out.println("<<<<<<<<<< DISTANCE MATRIX >>>>>>>>>>");
        System.out.println("*************************************\n");

        // prints the column header
        for (int column = 0; column<numCities; column++) {
            System.out.print("\t\t" + (column+1));
        }
        System.out.println("\n");
        
        for (int row = 0; row < numCities; row++) {
        	// prints the row headers
            System.out.print((row+1) + "\t");
            
            for (int column = 0; column < numCities; column++) {
            	
            	// if #row is not equal to #column, prints the distance between corresponding cities
            	// the double distances are rounded off to 6 decimal places
            	// else -> the main diagonal is 0, it is converted to int for easy interpretation
                if (row != column) {
                    System.out.printf("%.6f\t", distanceMatrix[row][column]);
                } else {
                    System.out.printf("%d\t\t", (int)distanceMatrix[row][column]);
                }
            }
            // line between each row
            System.out.println("\n");
        }
        System.out.println();
    }
    
    public static void setAllToUnvisited() {
    	// sets all elements of visited array to false
        for (int city = 0; city < numCities; city++) {
            visited[city] = false;
        }
    }
    
    public static void printVisitedList() {
    	// prints the visited array
        for (int city = 0; city < numCities; city++) {
            System.out.println(city + " " + visited[city]);
        }
    }
    
    public static void markVisited(int city) {
    	// markes "city" as visited in the "visited" array
        visited[city] = true;
    }
    
    public static void twoOptOptimizer() {
        ArrayList<Integer> currentRoute = new ArrayList<Integer>();
        ArrayList<Integer> newRoute = new ArrayList<Integer>();
        
        double newDistance;
        int improvement = INITIAL_IMPROVEMENT;
        
        // runs the nearest neighbor algorithm
        nearestNeighbor();
        
        // sets the "currentRoute" as the current "optimalRoute" (nearest neighbor route)
        currentRoute = optimalRoute;
        
        // the two opt optimizer algorithm, runs until no more improvement is made
    	while (improvement < IMPROVEMENT_LIMIT) {
        	// gets the distance of the "currentRoute"
            double currentBestDistance = calculateTotalDistance(currentRoute);
            /*  
             	the FIRST for loop starts from the second city in the list  
             	because the TSP requires the tour to start from a certain city 
             	and return to that certain city
             	otherwise it will mess up the order
             	runs up to the city before the last city
            */
            
            for (int cityA = 1; cityA < numCities - 1; cityA++) {
                for (int cityB = cityA + 1; cityB < numCities; cityB++) {
                	// gets the route and distance after the swap
                    newRoute = twoOptSwap(currentRoute, cityA, cityB);
                    newDistance = calculateTotalDistance(newRoute);
                	
                    // triggers if the new distance is better than the current best distance
                    if (newDistance < currentBestDistance) {
                    	// "improvement" is reset, "newRoute" becomes the "currentRoute" and "currentBestDistance" becomes "newDistance"
                        improvement = INITIAL_IMPROVEMENT;
                        currentRoute = newRoute;
                        currentBestDistance = newDistance;
                    }
                }
            }
            improvement++;
        }
    	
    	// set the current route from the two opt as the optimal route
    	optimalRoute = currentRoute;
        
        // calculates the total distance of the optimized route and prints the route and route length
        System.out.println("\n\n\n********************************************");
        System.out.println("<<<<<<<<<<<<<<< TWO OPT ROUTE >>>>>>>>>>>>>>");
        System.out.println("********************************************");
        
        double minimumDistance = calculateTotalDistance(currentRoute);
        printRoute(currentRoute);
        System.out.println("\n\nTWO OPT ROUTE LENGTH: " + minimumDistance);
        
    }
    
    public static ArrayList<Integer> twoOptSwap(ArrayList<Integer> route, int swapStart, int swapEnd) {
        
        ArrayList<Integer> newRoute = new ArrayList<Integer>();

        // adds all the cities up to the "swapStart" index to the newRoute
        for (int city = 0; city <= swapStart - 1; city++) {
            newRoute.add(route.get(city));
        }
        
        // adds the cities from "swapStart" to "swapEnd" in reverse order to the "newRoute"
        for (int city = swapEnd; city >= swapStart; city--) {
            newRoute.add(route.get(city));
        }
        
        // adds the remaining cities after "swapEnd" and adds them to the "newRoute"
        for (int city = swapEnd + 1; city < route.size(); city++) {
            newRoute.add(route.get(city));
        }

        // returns the "newRoute"
        return newRoute;
    }
    
    public static void nearestNeighbor() {
        double minimumDistance = MINIMUM_DISTANCE;
        int currentCity = FIRST_CITY;
        double minTracker = MAXIMUM_DISTANCE;
        
        // creates a route queue of MyQueue type
        MyQueue route = new MyQueue();

        // enqueues the first city to the route queue and sets it to visited
        route.enqueue(currentCity);
        visited[currentCity] = true;
        
        // marks the nextCity as not visited
        int nextCity = NOT_VISITED;
        
        // runs while the route queue is not empty 
        while (route.getSize() < numCities) {
        	
        	// loops from "city" to the total number of cities
            for (int city = 0; city < numCities; city++) {
            	// condition is met is "city" is unvisited and city is not the current city
            	if ((!visited[city]) && (city != currentCity)) {
                    
            		// condition is met if the distance from "currentCity" to "city" is less than "minTracker"
            		if (distanceMatrix[currentCity][city] <= minTracker) {
                        // "minTracker" is updated with the distance of "currentCity" to "city"
            			minTracker = distanceMatrix[currentCity][city];
                        // "nextCity" becomes "city"
            			nextCity = city;
                    }
                } 
            }
            
            // "nextCity" is enqueued in the route, gets marked as visited, "currentCity" becomes "nextCity"
            route.enqueue(nextCity);
            visited[nextCity] = true;
            currentCity = nextCity;
            // "minTracker" gets reinitialized with 999999
            minTracker = MAXIMUM_DISTANCE;
        }
   
        // prints the route and distance obtained from the nearest neighbor algorithm
        System.out.println("\n********************************************");
        System.out.println("<<<<<<<<<< NEAREST NEIGHBOR ROUTE >>>>>>>>>>");
        System.out.println("********************************************");

        // stores the first element in the current route queue
        int firstCityInRoute = route.peek();
        
        // calculates the total distance of the queue
        // runs while the route queue is not empty
        while (!route.isEmpty()) {
        	// takes the first city in the queue and adds it to the optimal route
            int previousCity = route.peek();
            optimalRoute.add(previousCity);
            
            // gets the city index and adds to make it readable (city 1 is indexed 0 in the array)
            System.out.print((previousCity + 1) + " -> ");
            // removes the front city of the route queue
            route.dequeue();
            
            // gets the distance from the "previousCity" to "next" adds it to the "minimumDistance"
            if (!route.isEmpty()) {
                int next = route.peek();
                minimumDistance += distanceMatrix[previousCity][next];
            } else {
            	// if the route queue is empty (for the last city in the route queue), "next" becomes the first city in the (completes the cycle)
                int next = firstCityInRoute;
                minimumDistance += distanceMatrix[previousCity][next];
            }
        }
        
        // adds the first city to the front of the route queue to complete the cycle
        optimalRoute.add(firstCityInRoute);
        // 1 is added to make it interpretable
        System.out.println(firstCityInRoute+1);
        // prints the total route length
        System.out.println("\nNEAREST NEIGHBOUR ROUTE LENGTH: " + minimumDistance);
        
        
        // returns to the twoOptOptimizer()
    }

    public static double calculateTotalDistance(ArrayList<Integer> route) {
        
    	// route contains the complete tour (first city is repeated as the last city)
        Object[] routeArray = route.toArray();
        
        double totalDistance = MINIMUM_DISTANCE;
               
        // gets the distance between two cities in the route and adds them in each iteration
        for (int city = 0; city < (routeArray.length - 1); city++) {
            totalDistance += distanceMatrix[(int)routeArray[city]][(int)routeArray[city + 1]];
        }
        
        // returns the total distance
        return totalDistance;
    }
    
    public static void printRoute(ArrayList<Integer> route) {
    	// converts the route ArrayList to an Object array
    	// route contains the complete tour (first city is repeated as the last city)
        Object[] routeArray = route.toArray();

    	// converts the object to int and prints all the cities
        for (int city = 0; city < routeArray.length; city++) {
        	// converts the object to int and prints all the cities, prints " -> " at the city until the last city
            System.out.print(((int)routeArray[city] + 1) + (city != (routeArray.length - 1) ? " -> " : ""));
        }
    }
    
    public static void printExecutionTime() {
    	// prints the execution time of the program in seconds
    	endTime = System.nanoTime();
        totalTime = endTime - startTime;
        double seconds = (double)totalTime / 1_000_000_000.0;
        
        System.out.println("\n\n***********************************************");
        System.out.println("<<<<< EXECUTION TIME: " + seconds + " SECONDS >>>>>");
        System.out.println("***********************************************\n");

    	// to get the time in nanoseconds, print the "totalTime" instead of the "seconds"
//        System.out.println("\n\n*************************************************");
//        System.out.println("<<<<< EXECUTION TIME: " + totalTime + " NANO SECONDS >>>>>");
//        System.out.println("*************************************************\n");
        
    }
    
    public static void dynamicProtrammingAlgorithm() {
    	// initializes/resets the optimal route
        optimalRoute = new ArrayList<Integer>();
        
        // runs the dynamic programming optimizer algorithm
        runDynamicProrgammingOptimizer();

        System.out.println("\n\n\n********************************************");
        System.out.println("<<<<<<<<< DYNAMIC PROGRAMMING ROUTE >>>>>>>>");
        System.out.println("********************************************");
        
        // prints the optimal route
        printRoute(optimalRoute);
        
        // prints the optimal distance
        System.out.print("\n\nDYNAMIC PROGRAMMING ROUTE LENGTH: " + optimalDistance);
    }

    
    public static void runDynamicProrgammingOptimizer() {
    	
    	// left shifts the binary value of 1 by numCities and decreases 1
    	// i.e numCities = 11; 1 << 11 equals 100000000000, which is equivalent to 1024 (2^11)
    	final int END_STATE = (1 << numCities) - 1;
    	
    	// creates a numCities x (2^numCities) costMatrix
    	Double[][] costMatrix = new Double[numCities][1 << numCities];

    	// adds all outgoing edges from the second city to the costMatrix
    	for (int city = 1; city < numCities; city++) {
    		costMatrix[city][(1 << FIRST_CITY) | (1 << city)] = distanceMatrix[FIRST_CITY][city];
    	}

    	for (int cityA = 3; cityA <= numCities; cityA++) {
    		for (int subset : combinations(cityA, numCities)) {
    			if (contains(FIRST_CITY, subset)) {
    				for (int nextCity = 0; nextCity < numCities; nextCity++) {
  	  	        		if (nextCity == FIRST_CITY || !contains(nextCity, subset)) continue;
          				int subsetWithoutNext = subset ^ (1 << nextCity);
  	          			double minimumDistance = Double.POSITIVE_INFINITY;
          				for (int city = 0; city < numCities; city++) {
  	          				if (city == FIRST_CITY || city == nextCity || !contains(city, subset)) {
  	          					continue;
  	          				}
  	          				double newDistance = costMatrix[city][subsetWithoutNext] + distanceMatrix[city][nextCity];
  	      					if (newDistance < minimumDistance) {
  	      						minimumDistance = newDistance;
  	      					}
  	  		          	}
  	  		          	costMatrix[nextCity][subset] = minimumDistance;
  	  	        	}
    			}
  	        
    		}
    	}

    	// connects the optimalRoute back to starting city and minimizes total distance
		for (int city = 1; city < numCities; city++) {
    		double currentDistance = costMatrix[city][END_STATE] + distanceMatrix[city][FIRST_CITY];
    		if (currentDistance < optimalDistance) {
    			optimalDistance = currentDistance;
    		}
		}

		int lastCity = FIRST_CITY;
		int state = END_STATE;
		optimalRoute.add(FIRST_CITY);

		// reconstructs the route using the costMatrix
		for (int cityA = 1; cityA < numCities; cityA++) {
        
			int currentCity = -1;
			for (int cityB = 0; cityB < numCities; cityB++) {
				if (cityB == FIRST_CITY || !contains(cityB, state)) {
					continue;
				}
				if (currentCity == -1) {
					currentCity = cityB;
				}
				double prevDist = costMatrix[currentCity][state] + distanceMatrix[currentCity][lastCity];
				double newDist  = costMatrix[cityB][state] + distanceMatrix[cityB][lastCity];
				if (newDist < prevDist) {
					currentCity = cityB;
				}
			}

			optimalRoute.add(currentCity);
			state = state ^ (1 << currentCity);
			lastCity = currentCity;
		}

		// adds the first city to complete the tour
		optimalRoute.add(FIRST_CITY);
      
		// reverses the route
		Collections.reverse(optimalRoute);

	}

    private static boolean contains(int element, int subset) {
    	// returns true if the city is present in the subset
    	return ((1 << element) & subset) != 0;
    }

    // generates all bit sets of size where number of cities bits are set to one. subsets are returned as a list of integer masks
    public static ArrayList<Integer> combinations(int city, int size) {
    	ArrayList<Integer> subsets = new ArrayList<>();
    	combinations(0, 0, city, size, subsets);
    	return subsets;
    }

    private static void combinations(int set, int position, int currentCity, int size, List<Integer> subsets) {
      
    	// return early if there are more cities left to select than what is available
    	int elementsLeftToPick = size - position;
    	
    	if (elementsLeftToPick < currentCity) {
    		return;
    	}

    	// 'currentCity' elements are selected, a valid subset is found
    	if (currentCity == 0) {
    		// adds the set to the subset list
    		subsets.add(set);
    	} else {
    		for (int city = position; city < size; city++) {
    			// tries to insert this city to the set
    			set ^= (1 << city);

    			combinations(set, city + 1, currentCity - 1, size, subsets);

    			// Backtracks and tries the instance where this city was not included
    			set ^= (1 << city);
    		}
		}
	}
    
}
    
