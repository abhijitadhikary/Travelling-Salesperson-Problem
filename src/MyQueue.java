public class MyQueue {

	// this queue is only meant for int type elements
    private int[] array;
    private int size;
    
    // constructor for the queue, sets the size = 0
    public MyQueue() {
        this.size = 0;
    }

    public boolean enqueue(int element) {
    	// if the queue is not initialized, initializes it with a 1x1 int array
    	if (this.size == 0) {
            array = new int[1];
        }
        
    	// increases the array size by 1
        increaseArraySize();
        
        // stores the element at the last position of the array
        this.array[size-1] = element;
        
        // returns true if insertion is successful
        return true;
    }
    
    public boolean dequeue() {
        if (!isEmpty()) {
        	// removes the first element of the array
            decreaseArraySize();
            return true;
        } else {
            return false;
        }
    }
    
    public int peek() {
    	// returns the first element of the array
        if (!isEmpty()) {
            return this.array[0];
        } else {
            return -1;
        }
    }
    
    public boolean isEmpty() {
    	// returns true if the array is empty
        return this.size <= 0; 
    }
    
    public void increaseArraySize() {
    	// increases the array size by one and copies all the elements from the old array to the new array, then updates the current array and the current size
        int newSize = this.size + 1;
        int[] newArray = new int[newSize];
        
        System.arraycopy(this.array, 0, newArray, 0, this.array.length);
        
        this.array = newArray;
        this.size = newSize;
    }
    
    public void decreaseArraySize() {
    	// creates a new array of size-1 and copies all the elements(except the first one) of the existing array to the new array, then assigns the new array as the current array
        int newSize = this.size - 1;
        int[] newArray = new int[newSize];
        
        System.arraycopy(this.array, 1, newArray, 0, newSize);
        
        this.array = newArray;
        this.size = newSize;
    }

    
    public void printQueue() {
    	// prints the whole array according to insertion order
        for (int i: this.array) {
            System.out.print(i + " ");
        }
    }
    
    public int getSize() {
    	// returns the size of the array if it is not empty
        if (!isEmpty()) {
            return this.size;
        } else {
            return -1;
        }
    }
    
    public void shuffle() {
    	if (!isEmpty()) {
    		int temp = array[size-1];
    		
    		for (int i = 1; i<array.length; i++) {
    			array[i] = array[i-1];
    		}
    		
    		array[0] = temp;
    		
    	}
    }
    
    
}