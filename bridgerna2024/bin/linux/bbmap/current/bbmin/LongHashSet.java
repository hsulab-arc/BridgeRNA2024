package bbmin;

import java.util.Arrays;
import java.util.Random;

/**
 * @author Brian Bushnell
 * @date July 6, 2016
 *
 */
public final class LongHashSet{
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	public LongHashSet(){
		this(256);
	}
	
	public LongHashSet(int initialSize){
		this(initialSize, 0.7f);
	}
	
	public LongHashSet(int initialSize, float loadFactor_){
		invalid=randy.nextLong()|MINMASK;
		assert(invalid<0);
		assert(initialSize>0) : "Attempting to initialize a "+getClass().getSimpleName()+" of size<1.";
		assert(loadFactor_>0 && loadFactor_<1) : "Attempting to initialize a "+getClass().getSimpleName()+" with invalid load factor: "+loadFactor_;
		loadFactor=mid(0.25f, loadFactor_, 0.90f);
		resize(initialSize);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Public Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	public void clear(){
		if(size<1){return;}
		Arrays.fill(array, invalid);
		size=0;
	}
	
	public boolean contains(long value){
		return value==invalid ? false : findCell(value)>=0;
	}
	
	/**
	 * Add this value to the set.
	 * @param value
	 * @return true if the value was added, false if it was already contained.
	 */
	public boolean add(long value){
		if(value==invalid){resetInvalid();}
		int cell=findCellOrEmpty(value);
		if(array[cell]==invalid){
			array[cell]=value;
			size++;
			if(size>sizeLimit){resize();}
			return true;
		}
		assert(array[cell]==value);
		return false;
	}
	
	/**
	 * Remove this value from the set.
	 * @param value
	 * @return true if the value was removed, false if it was not present.
	 */
	public boolean remove(long value){
		if(value==invalid){return false;}
		final int cell=findCell(value);
		if(cell<0){return false;}
		assert(array[cell]==value);
		array[cell]=invalid;
		size--;
		
		rehashFrom(cell);
		return true;
	}
	
	public int size(){return size;}
	
	public boolean isEmpty(){return size==0;}
	
	/*--------------------------------------------------------------*/
	/*----------------        String Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	public String toString(){
		return toStringListView();
	}
	
	public String toStringSetView(){
		StringBuilder sb=new StringBuilder();
		sb.append('[');
		String comma="";
		for(int i=0; i<array.length; i++){
			if(array[i]!=invalid){
				sb.append(comma+"("+i+", "+array[i]+")");
				comma=", ";
			}
		}
		sb.append(']');
		return sb.toString();
	}
	
	public String toStringListView(){
		StringBuilder sb=new StringBuilder();
		sb.append('[');
		String comma="";
		for(int i=0; i<array.length; i++){
			if(array[i]!=invalid){
				sb.append(comma+array[i]);
				comma=", ";
			}
		}
		sb.append(']');
		return sb.toString();
	}
	
	public long[] toArray(){
		long[] x=new long[array.length];
		int i=0;
		for(long v : array){
			x[i]=v;
			i++;
		}
		return x;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Private Methods       ----------------*/
	/*--------------------------------------------------------------*/
	
	public boolean verify(){
		int numValues=0;
		int numFound=0;
		for(int i=0; i<array.length; i++){
			final long value=array[i];
			if(value!=invalid){
				numValues++;
				final int cell=findCell(value);
				if(i==cell){
					numFound++;
				}else{
					return false;
				}
			}
		}
		return numValues==numFound && numValues==size;
	}
	
	private void rehashFrom(int initial){
		if(size<1){return;}
		final int limit=array.length;
		for(int cell=initial+1; cell<limit; cell++){
			final long x=array[cell];
			if(x==invalid){return;}
			rehashCell(cell);
		}
		for(int cell=0; cell<initial; cell++){
			final long x=array[cell];
			if(x==invalid){return;}
			rehashCell(cell);
		}
	}
	
	private boolean rehashCell(final int cell){
		final long value=array[cell];
		assert(value!=invalid);
		if(value==invalid){resetInvalid();}
		final int dest=findCellOrEmpty(value);
		if(cell==dest){return false;}
		assert(array[dest]==invalid);
		array[cell]=invalid;
		array[dest]=value;
		return true;
	}
	
	private void resetInvalid(){
		final long old=invalid;
		long x=invalid;
		while(x==old || contains(x)){x=randy.nextLong()|MINMASK;}
		assert(x<0);
		invalid=x;
		for(int i=0; i<array.length; i++){
			if(array[i]==old){array[i]=invalid;}
		}
	}
	
	private int findCell(final long value){
		if(value==invalid){return -1;}
		
		final int limit=array.length, initial=(int)((value&MASK)%modulus);
		for(int cell=initial; cell<limit; cell++){
			final long x=array[cell];
			if(x==value){return cell;}
			if(x==invalid){return -1;}
		}
		for(int cell=0; cell<initial; cell++){
			final long x=array[cell];
			if(x==value){return cell;}
			if(x==invalid){return -1;}
		}
		return -1;
	}
	
	private int findCellOrEmpty(final long value){
		assert(value!=invalid) : "Collision - this should have been intercepted.";
		
		final int limit=array.length, initial=(int)((value&MASK)%modulus);
		for(int cell=initial; cell<limit; cell++){
			final long x=array[cell];
			if(x==value || x==invalid){return cell;}
		}
		for(int cell=0; cell<initial; cell++){
			final long x=array[cell];
			if(x==value || x==invalid){return cell;}
		}
		throw new RuntimeException("No empty cells - size="+size+", limit="+limit);
	}
	
	public final void resizeDestructive(int newSize){
		size=0;
		sizeLimit=0;
		array=null;
		resize(newSize);
	}
	
	private final void resize(){
		assert(size>=sizeLimit);
		resize(array.length*2L+1);
	}
	
	private final void resize(final long size2){
		assert(size2>size) : size+", "+size2;
		
		//This is supposed to be a prime but the primes code is ripped out in this version.
		//Any odd number is fine in most cases.
		long newPrime=size2|1;
		if(newPrime+extra>Integer.MAX_VALUE){
			newPrime=(Integer.MAX_VALUE-extra-2)|1;
		}
		assert(newPrime>modulus) : "Overflow: "+size+", "+size2+", "+modulus+", "+newPrime;
		modulus=(int)newPrime;
		
		final int size3=(int)(newPrime+extra);
		sizeLimit=(int)(modulus*loadFactor);
		final long[] old=array;
		array=new long[size3];
		Arrays.fill(array, invalid);
		
//		System.err.println("Resizing "+(old==null ? "null" : ""+old.length)+" to "+size3);
		
		if(size<1){return;}
		
		size=0;
		for(long value : old){
			if(value!=invalid){
				add(value);
			}
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------      Stuff From BBTools      ----------------*/
	/*--------------------------------------------------------------*/
	
	public static final float mid(float x, float y, float z){
		return x<y ? (x<z ? min(y, z) : x) : (y<z ? min(x, z) : y);
	}
	public static final float min(float x, float y){return x<y ? x : y;}
	public static final float max(float x, float y){return x>y ? x : y;}
	
	/** Number of values that can be held without resizing */
	public int capacity(){
		return sizeLimit;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	private long[] array;
	private int size=0;
	/** Value for empty cells */
	private long invalid;
	private int modulus;
	private int sizeLimit;
	private final float loadFactor;
	
	private static final Random randy=new Random(1);
	private static final long MASK=Long.MAX_VALUE;
	private static final long MINMASK=Long.MIN_VALUE;
	
	private static final int extra=10;
	
}
