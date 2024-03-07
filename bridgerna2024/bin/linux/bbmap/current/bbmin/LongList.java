package bbmin;

import java.util.Arrays;

public final class LongList{
	
	public LongList(){this(256);}
	
	public LongList(int initial){
		assert(initial>0);
		array=allocLong1D(initial);
	}
	
	public void clear(){
		size=0;
	}
	
	public final void set(int loc, long value){
		if(loc>=array.length){
			resize(loc*2L+1);
		}
		array[loc]=value;
		size=max(size, loc+1);
	}
	
	public final void setLast(long value){
		assert(size>0);
		array[size-1]=value;
	}
	
	public final void increment(int loc, long value){
		if(loc>=array.length){
			resize(loc*2L+1);
		}
		array[loc]+=value;
		size=max(size, loc+1);
	}
	
	public final void increment(int loc){
		increment(loc, 1);
	}
	
	public final void incrementBy(LongList b){
		for(int i=b.size-1; i>=0; i--){
			increment(i, b.get(i));
		}
	}
	
	public final void incrementBy(long[] b){
		for(int i=b.length-1; i>=0; i--){
			increment(i, b[i]);
		}
	}
	
	public final void append(LongList b){
		for(int i=0; i<b.size; i++){
			add(b.get(i));
		}
	}
	
	public final void append(long[] b){
		for(int i=0; i<b.length; i++){
			add(b[i]);
		}
	}
	
	public final long get(int loc){
		return(loc>=size ? 0 : array[loc]);
	}
	
	public final void add(long x){
		if(size>=array.length){
			resize(size*2L+1);
		}
		array[size]=x;
		size++;
	}
	
	private final void resize(final long size2){
		assert(size2>size) : size+", "+size2;
		final int size3=(int)min(MAX_ARRAY_LEN, size2);
		assert(size3>size) : "Overflow: "+size+", "+size2+" -> "+size3;
		array=copyOf(array, size3);
	}
	
	public final void shrink(){
		if(size==array.length){return;}
		array=copyOf(array, size);
	}
	
	public final double stdev(){
		if(size<2){return 0;}
		double sum=sum();
		double avg=sum/size;
		double sumdev2=0;
		for(int i=0; i<size; i++){
			long x=array[i];
			double dev=avg-x;
			sumdev2+=(dev*dev);
		}
		return Math.sqrt(sumdev2/size);
	}
	
	public final double avgDif(final double x){
		double sum=0;
		for(int i=0; i<size; i++){
			sum+=absdif(x, array[i]);
		}
		return sum/(max(1, size));
	}
	
	public final double rmsDif(final double x){
		double sum=0;
		for(int i=0; i<size; i++){
			double dif=absdif(x, array[i]);
			sum+=dif*dif;
		}
		return Math.sqrt(sum/(max(1, size)));
	}
	
	public final long sumLong(){
		long sum=0;
		for(int i=0; i<size; i++){
			sum+=array[i];
		}
		return sum;
	}
	
	public final double sum(){
		double sum=0;
		for(int i=0; i<size; i++){
			sum+=array[i];
		}
		return sum;
	}
	
	public final double mean(){
		return size<1 ? 0 : sum()/size;
	}
	
	//Ignores elements below 1
	public final double harmonicMean(){
		double sum=0;
		int count=0;
		for(int i=0; i<size; i++){
			if(array[i]>0){
				sum+=1.0/array[i];
				count++;
			}
		}
		double avg=sum/max(1, count);
		return 1.0/avg;
	}
	
	//Ignores elements below 1
	public final double geometricMean(){
		double sum=0;
		int count=0;
		for(int i=0; i<size; i++){
			if(array[i]>0){
				sum+=Math.log(array[i]);
				count++;
			}
		}
		double avg=sum/max(1, count);
		return Math.exp(avg);
	}
	
	/** Assumes list is sorted */
	public final double medianWeightedAverage(){
		if(size<1){return 0;}
		int half=size/2;
		long count=0;
		double sum=0;
		for(int i=0, j=size-1; i<half; i++, j--){
			int mult=i+1;
			double incr=(array[i]+array[j])*mult;
			sum+=incr;
			count+=2*mult;
		}
		if((size&1)==1){//odd length
			int mult=half+1;
			double incr=(array[half])*mult;
			sum+=incr;
			count+=2*mult;
		}
		return sum/count;
	}
	
	/** Assumes list is sorted */
	public final long median(){
		if(size<1){return 0;}
		int idx=percentileIndex(0.5);
		return array[idx];
	}
	
	/** Allows unsorted list */
	public final long min(){
		if(size<1){return 0;}
		long x=array[0];
		for(int i=1; i<size; i++){
			x=min(x, array[i]);
		}
		return x;
	}
	
	/** Allows unsorted list */
	public final long max(){
		if(size<1){return 0;}
		long x=array[0];
		for(int i=1; i<size; i++){
			x=max(x, array[i]);
		}
		return x;
	}
	
	/** Assumes list is sorted */
	public final long mode(){
		if(size<1){return 0;}
		assert(sorted());
		int streak=1, bestStreak=0;
		long prev=array[0];
		long best=prev;
		for(int i=0; i<size; i++){
			long x=array[i];
			if(x==prev){streak++;}
			else{
				if(streak>bestStreak){
					bestStreak=streak;
					best=prev;
				}
				streak=1;
				prev=x;
			}
		}
		if(streak>bestStreak){
			bestStreak=streak;
			best=prev;
		}
		return best;
	}
	
	public long percentile(double fraction){
		if(size<1){return 0;}
		int idx=percentileIndex(fraction);
		return array[idx];
	}
	
	public int percentileIndex(double fraction){
		if(size<2){return size-1;}
		assert(sorted());
		double target=(sum()*fraction);
		double sum=0;
		for(int i=0; i<size; i++){
			sum+=array[i];
			if(sum>=target){
				return i;
			}
		}
		return size-1;
	}
	
	public final void shrinkToUnique(){
		condense();
		shrink();
	}
	
	//In-place.
	//Assumes sorted.
	public final void condense(){
		if(size<=1){return;}
		
		int i=0, j=1;
		for(; j<size && array[i]<array[j]; i++, j++){}//skip while strictly ascending 
		
		int dupes=0;
		for(; j<size; j++){//This only enters at the first nonascending pair
			long a=array[i], b=array[j];
			assert(a<=b) : "Unsorted: "+i+", "+j+", "+a+", "+b;
			if(b>a){
				i++;
				array[i]=b;
			}else{
				//do nothing
				dupes++;
				assert(a==b);
			}
		}
		assert(dupes==(size-(i+1)));
		assert(size>=(i+1));
		size=i+1;
	}
	
	@Override
	public String toString(){
		return toStringListView();
	}
	
	public String toStringSetView(){
		StringBuilder sb=new StringBuilder();
		sb.append('[');
		String comma="";
		for(int i=0; i<size; i++){
			if(array[i]!=0){
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
		for(int i=0; i<size; i++){
				sb.append(comma+array[i]);
				comma=", ";
		}
		sb.append(']');
		return sb.toString();
	}
	
	public long[] toArray(){
		long[] x=allocLong1D(size);
		for(int i=0; i<x.length; i++){
			x[i]=array[i];
		}
		return x;
	}
	
	public void sort() {
		if(size>1){Arrays.sort(array, 0, size);}
	}
	
	public void sortSerial() {
		if(size>1){Arrays.sort(array, 0, size);}
	}
	
	public void reverse() {
		if(size>1){reverseInPlace(array, 0, size);}
	}
	
	public boolean sorted(){
		for(int i=1; i<size; i++){
			if(array[i]<array[i-1]){return false;}
		}
		return true;
	}
	
	public int size() {
		return size;
	}
	
	public int capacity() {
		return array.length;
	}
	
	public int freeSpace() {
		return array.length-size;
	}
	
	private static void reverseInPlace(final long[] array, final int from, final int to){
		if(array==null){return;}
		final int len=to-from;
		final int max=from+len/2, last=to-1;
		for(int i=from; i<max; i++){
			long temp=array[i];
			array[i]=array[last-i];
			array[last-i]=temp;
		}
	}
	
	private static final long min(long x, long y){return x<y ? x : y;}
	private static final long max(long x, long y){return x>y ? x : y;}
	
	private static final int min(int x, int y){return x<y ? x : y;}
	private static final int max(int x, int y){return x>y ? x : y;}
	
	private static double absdif(double a, double b) {return a>b ? a-b : b-a;}
	
	private static final long[] allocLong1D(int x){return new long[x];}
	private static long[] copyOf(long[] buffer, long newLength) {
		final int len=buffer.length;
		final int len2=(int)min(newLength, MAX_ARRAY_LEN);
		if(newLength>len2 && len2<=len){
			throw new RuntimeException("Tried to create an array above length limit: "+len+"," +newLength);
		}
		long[] copy=null;
		try {
			copy=Arrays.copyOf(buffer, len2);
		} catch (OutOfMemoryError e) {
			throw e;
		}
		return copy;
	}
	
	private static final long MAX_ARRAY_LEN=Integer.MAX_VALUE-20;
	
	public long[] array;
	/** Highest occupied index plus 1, i.e., lowest unoccupied index */
	public int size=0;
	
}
