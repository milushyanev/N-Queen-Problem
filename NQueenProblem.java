/*
Milush Yanev
10.2.2019
CS 4200
N-Queen Problem
*/

package nqueenproblem;
import java.util.Arrays;
import java.util.Comparator;
import java.util.Scanner;

/**
 *
 * @author milus
 */

class GeneticAlgorithm {


    public int[] solve(int n, int populationSize, double mProb, int gen) {

        // each one should get a mate.
        populationSize = populationSize - (populationSize % 2); 

        //generate current population (queens,selected size)
        int[][] population = genPop(n, populationSize);

        //get the maximum fitness that the queens can get
        int mF = get_mF(n);

        //loop until number of generations is reached
        for (int i = 0; i < gen; i++) {
            
            //get the population and handle crossover
            population = getPopulation(population);

            population = CrossoverH(population, n);

            for (int j = 0; j < populationSize; j++) {

                //get the fitness by setting the current element to max fitness
                if (getFitness(population[j]) == mF)
                    return population[j];

                //mutate by getting population index and max probability
                population[j] = mutation(population[j], mProb);
                
                //if condition is met return current index
                if (getFitness(population[j]) == mF)
                    return population[j];

            }

        }
        return null;
    }

    //handle crossover
    private int[][] CrossoverH(int[][] population, int n) {
        for (int i = 0; i < population.length; i += 2) {
            
            //set the position
            int crossoverPos = (int) (Math.random() * n);

            //loop and swap until crossover is handled
            for (int j = 0; j < crossoverPos; j++) {
                int tmp = population[i][j];
                population[i][j] = population[i+1][j];
                population[i+1][j] = tmp;
            }

        }
        //return population after crossover is handled
        return population;
    }
    
    //get the population by using lambda and Comparator function
    private int[][] getPopulation(int[][] population) {
        //compare and sort population
        Arrays.sort(population, Comparator.comparingInt(this::getFitness));

        return population;
    }
    
    //try to mutate
    private int[] mutation(int[] res, double mutationProbability) {
        //if the probability is satisfied mutate
        if (satisfyProb(mutationProbability))
            res[(int)(Math.random()*res.length)] = (int)(Math.random()*res.length);
        
        //return the result
        return res;
    }

    //set probability condition
    private boolean satisfyProb(double prob) {
        return prob >= Math.random();
    }
    
    //get Fitness
    private int getFitness(int[] res) {
        //fitness=max fitness-heuristic cost
        return get_mF(res.length) - Solve.getHeuristicCost(res);
    }
    
    //get max fitness by using the formula 
    private int get_mF(int n) {
        return n*(n-1)/2;
    }

    //generate the solution by setting to random state
    private int[] genSolution(int n) {
        return Solve.setState(n);
    }

    //get the generation probability
    private int[][] genPop(int n, int nPop) {
        int[][] population = new int[nPop][];
        //loop until the current population is reached
        for (int i = 0; i < nPop; i++)
            
            //set proper index for population and generated solution
            population[i] = genSolution(n);

        return population;
    }
}

  //class Solver to handle both algorithms
  class Solve {
    // Generate state that all queens have row # 0
    public static int[] oneDqueens(int n) {

        return new int[n];
    }
    // Randomizes state
    public static int[] randomState(int[] res) {

        for (int i = 0; i < res.length; i++)
            res[i] = (int) (Math.random() * res.length);

        return res;
    }
    // generates random initial state
    public static int[] setState(int n) {
        return randomState(oneDqueens(n));
    }
    // Returns heuristic cost
    public static int getHeuristicCost(int[] res) {
        int h = 0;
        // increment cost if two queens are in same row or in same diagonal.
        for (int i = 0; i < res.length; i++)
            for (int j = i + 1; j < res.length; j++)
                if (res[i] == res[j] || Math.abs(res[i] - res[j]) == j - i)
                    h += 1;
        return h;
    }
 }

class SAnnealing {

    //set conditions for Simulated Annealing solving - board size,maxIterations,temperature,cooling Factor
    public int[] solve(int boardSize, int mIter, double temp, double cFactor) {
        int[] result = Solve.setState(boardSize);

        //get the cost to beat by getting heuristic
        int costToBeat = Solve.getHeuristicCost(result);

        // terminate when it reaches max num of iterations or problem is solved.
        for (int x = 0; x < mIter && costToBeat > 0; x++) {
            result = move(result, costToBeat, temp);
            costToBeat = Solve.getHeuristicCost(result);
            //current temp = temp*cooling factor , 0.01
            temp = Math.max(temp * cFactor, 0.01);
        }
        // return solution if solved
        return costToBeat == 0 ? result : null; 
    }

    //check moves
    private int[] move(int res[], int costToBeat, double temp) {
        int n = res.length;

        while (true) {
            int i = (int) (Math.random() * n);
            int j = (int) (Math.random() * n);
            int tmpRow = res[i];
            res[i] = j;

            int cost = Solve.getHeuristicCost(res);
            if (cost < costToBeat)
                return res;
            
            //delta E = cost to beat minus current cost
            int dE = costToBeat - cost;
            double acceptProb = Math.min(1, Math.exp(dE / temp));

            if (Math.random() < acceptProb)
                return res;

            res[i] = tmpRow;
        }
    }

}

public class NQueenProblem {

//print the solution
public static void printSolution(int[] res,int n, int board[][]){
//message if result is returning null
if(res==null){
     System.out.print("No solution found ");
     }else 
{
      //loop until result index matches column/row index
      for (int i=0;i<n;i++){

          for (int j=0;j<n;j++){
              //if current index is empty return 0
              if(j!=res[i]){
              board[i][j]=0;
              //if current index is not empty return a queen or 1 
              }else
              board[i][j]=1;
          }
      }
 //for loop to print result in 2d array
 for (int i = 0; i < n; i++) { 
            for (int j = 0; j < n; j++) 
                System.out.print(" " + board[i][j] 
                                 + " "); 
            System.out.println(); 
        } 
}
}
//start N queen menu options
public static void startGame(){
		System.out.println("Select an option to start N Queen Proble ");
		System.out.println("1. Simulated Annealing");
		System.out.println("2. GeneticAlgorithm");
		System.out.println("3. EXIT");
		System.out.print("\n->");
}

    public static void main(String[] args) {
        
         int choice ;
         SAnnealing simulatedAnnealing = new SAnnealing();
         GeneticAlgorithm geneticAlgorithm = new GeneticAlgorithm();
         //number of queens
         int N=10;
         
         //number of Simulated Annealing Itterations
         int saIterations=50000;
         
         //set temperature
         double saTemperature=120;
         
         //set cooling factor
         double saCoolingFactor=.95;
         
         //Genetic Algorithm population
         int nPop=2000;
         
         //get the mutation probability for GA
         double mProb=0.2;
         
         //initiate the number of generations for GA
         int nGen=50000;
         
         //set board to 2d
         int board[][]=new int[N][N];
         
        //call menu
        startGame();
        try (Scanner sc = new Scanner(System.in)) {
            choice=sc.nextInt(3);
            if (choice==1) {
                long startTime = System.nanoTime();
                int[] res = simulatedAnnealing.solve(N, saIterations, saTemperature, saCoolingFactor);
                printSolution(res,N,board);
                long endTime = System.nanoTime();
                long timeElapsed = endTime - startTime;
                System.out.println(" ");
                System.out.println("Number of queens on board set to  "+N);
                System.out.println("Max Number of Iterations set to  "+saIterations);
                System.out.println("Execution time in milliseconds : " +
                        timeElapsed / 1000000);
                System.out.println(" ");
            }
            else if(choice==2) {
                
                long startTime = System.nanoTime();
                int[] res= geneticAlgorithm.solve(N,nPop,mProb,nGen);
                printSolution(res,N,board);
                long endTime = System.nanoTime();
                long timeElapsed = endTime - startTime;
                
                System.out.println(" ");
                System.out.println("Population set to " + nPop);
                System.out.println("Mutation probability is " + mProb);
                System.out.println("Number of generations is " + nGen);
                System.out.println("Execution time in milliseconds : " +
                        timeElapsed / 1000000);
                System.out.println(" ");
            }
            else if(choice==3) {
                System.exit(0);
            }
        }     
    }
    
}
