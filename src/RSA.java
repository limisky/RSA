import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.math.BigInteger;
import java.util.Random;

public class RSA {
	final static int L = 3;					//size of block - the encryption deals with L characters at a time
	final static int PQ_BITLENGTH = 32;		//size of p,q
	final static int KEY_BITLENGTH = 32;	//size of e
	
	static Random rnd = new Random();
	static BigInteger big_two = new BigInteger("2");
	
	public static void main(String[] args) {
		
		//System.out.println("test sqm: "+sqm(new BigInteger("3"),new BigInteger("5"),new BigInteger("7")).toString());
		//System.out.println("test generatePrime: "+generatePrime(PRIME_BITLENGTH));
		
		long startTime = System.nanoTime();
		
		BigInteger[] keyPair = generateKey();
		
		//System.out.println("public key: (" + keyPair[0] + "," + keyPair[1] + ")");
		//System.out.println("private key: (" + keyPair[2] + "," + keyPair[3] + ")");
		//System.out.println("test keyPair: "+sqm(sqm(big_two,keyPair[0],keyPair[1]),keyPair[2],keyPair[3]));
		
		BigInteger[] cipher = encrypt(readFile("rsa_group9_"+L+".plain"),L,keyPair[0],keyPair[1]);
		
		long endTime = System.nanoTime();
		System.out.println("Run time： "+(endTime-startTime)+"ns"); 
		
		writeFile(cipher,keyPair);
	}
	/*
	 * Write cipher into file "rsa_group9_L.crypto"
	 * Write private key into file "rsa_group9_L.key"
	 * Write public key into file "rsa_group9_L.pub"
	 */
	static void writeFile(BigInteger[] cipher, BigInteger[] keyPair)
	{
		try{
			//crypto
			File crypto = new File("rsa_group9_"+L+".crypto");
			crypto.createNewFile();
			FileWriter cryptoWriter = new FileWriter(crypto);
			for(int i=0;cipher[i]!=null;i++)
				cryptoWriter.write(cipher[i].toString()+"\n");
			cryptoWriter.close();
			
			//private key
			File key = new File("rsa_group9_"+L+".key");
			key.createNewFile();
			FileWriter keyWriter = new FileWriter(key);
			keyWriter.write(keyPair[2].toString()+"\n");
			keyWriter.write(keyPair[3].toString()+"\n");
			keyWriter.close();
			
			//public key
			File pub = new File("rsa_group9_"+L+".pub");
			pub.createNewFile();
			FileWriter pubWriter = new FileWriter(pub);
			pubWriter.write(keyPair[0].toString()+"\n");
			pubWriter.write(keyPair[1].toString()+"\n");
			pubWriter.close();
			
		}catch(Exception e){
			e.printStackTrace();
		}
	}
	//e(k) = x^e mod n where x represent L characters
	static BigInteger[] encrypt(int[] plain, int L, BigInteger key_e, BigInteger key_n)
	{
		BigInteger[] cipher= new BigInteger[100/L+1];
		int temp = 0;
		int index = 0;
		for(int i=0,l=L;plain[i]!=-1;i++)
		{
			temp = (temp<<8) + plain[i];
			if(l==1)
			{
				cipher[index] = sqm(new BigInteger(""+temp),key_e,key_n);//encrypt
				l=L;
				index++;
				temp = 0;
			}
			else
				l--;
		}
		return cipher;
	}
	/*
	 *  Read plain text into an integer array,
	 *  characters are encoded on 8 bits (typically Latin1 encoding)
	 */
	static int[] readFile(String fileName)
	{
		int[] plain = new int[102];
		for(int i=0;i<plain.length;i++)
			plain[i]=-1;
		int temp = 0;
		int i = 0;
		try {
			FileReader reader = new FileReader(fileName);
			while((temp=reader.read())!=-1)
			{
				plain[i] = temp;
				i++;
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
		//when the number of characters to encode is not a multiple of L
		//------------------------------------------------------------
		//----------------???WHAT TO DO HERE???-----------------------
		//------------------------------------------------------------
		while(i%L!=0)
		{
			plain[i]=0;
			i++;
		}
        return plain;
	}
	/*
	 * Generate key pair for RSA
	 * @return 
	 * 		result[0]&result[1] - public key
	 * 		result[2]&result[3] - private key
	 */
	static BigInteger[] generateKey()
	{
		BigInteger result[] = new BigInteger[4];
		BigInteger p = generatePrime(PQ_BITLENGTH);
		BigInteger q,phi,n,e,d;
		do{
			q = generatePrime(PQ_BITLENGTH);
		}while (q.equals(p));
		n = p.multiply(q);//n = p*q
		phi = (p.subtract(BigInteger.ONE)).multiply((q.subtract(BigInteger.ONE)));//phi = (p-1)*(q-1)
		do{
			e = generatePrime(KEY_BITLENGTH);
		}while(!gcd(e,phi).equals(BigInteger.ONE)||e.compareTo(phi)>=0);
		d = inv(e,phi);
		
		result[0]=e;
		result[1]=n;
		result[2]=d;
		result[3]=n;
		
		return result;
	}
	/*
	 * Use Extended Euclidean Algorithm to find y = inv(x,n) is such that (x*y=1) mod n
	 * ax+by=1 return x+b
	 */
	static BigInteger inv(BigInteger x, BigInteger n)
	{
		return ext_gcd(x,n)[0].add(n);
	}
	static BigInteger[] ext_gcd(BigInteger a, BigInteger b)
	{
		BigInteger [] tempxy = new BigInteger[2];
		BigInteger temp = null;
		if (b.equals(BigInteger.ZERO))
		{
			tempxy[0]=BigInteger.ONE;
			tempxy[1]=BigInteger.ZERO;
		}
		else
		{
			tempxy = ext_gcd(b,a.mod(b));
			temp = tempxy[0];
			tempxy[0] = tempxy[1];
			tempxy[1] = temp.subtract((tempxy[1].multiply(a.divide(b))));
		}
		return tempxy;		
	}
	
	//calculate the greatest common divisor of a and b
	static BigInteger gcd(BigInteger a, BigInteger b)
	{
		BigInteger temp;
		while(!b.equals(BigInteger.ZERO))
		{
			temp = b;
			b = a.mod(b);
			a = temp;
		}
		return a;
	}
	//generate a strong probable prime with the specified bitLength
	static BigInteger generatePrime(int bitLength){
		BigInteger x = new BigInteger(bitLength,rnd);//[0,2^bitLength-1]
		//change the most significant bit to 1 to make sure it meets the bitLength requirement
		x = x.or(BigInteger.ONE.shiftLeft(bitLength-1));//[2^(bitLength-1),2^bitLength-1]
		if(x.mod(big_two).equals(BigInteger.ZERO))//if n is even
			x = x.add(BigInteger.ONE);
		while(!miller_rabin(x,5))
			x = x.add(big_two);
		return x;
	}
	/**
	 * Use Miller-Rabin to test if n is a strong probable prime
	 * @param n - number to be tested from primality
	 * @param k - How many base a to be tested
	 * @return
	 * 		 true means n has 1-4^(-k) probability to be prime
	 * 		false means n is composite
	 */
	static boolean miller_rabin(BigInteger n, int k)
	{
		BigInteger a;
        for (;k>0;k--) {
        	//pick a random integer a in the range [1,N-1]
            do{
                a = new BigInteger(n.bitLength(), rnd);
            }while(a.equals(BigInteger.ZERO)||a.compareTo(n)>=0);
            if (!miller_rabin_test(a, n))
                return false;
        }
        return true;	
	}
	//Test if n is a strong probable prime to base a
    static boolean miller_rabin_test(BigInteger a, BigInteger n) {
    	BigInteger n_minus_one = n.subtract(BigInteger.ONE);
    	BigInteger d = n_minus_one;
		int s = d.getLowestSetBit();
		d = d.shiftRight(s);
		BigInteger x = sqm(a,d,n);
		if(x.equals(BigInteger.ONE)||x.equals(n_minus_one))
			return true;
		for(int i=0;i<s-1;i++)
		{
			x = x.multiply(x).mod(n);	
			if(x.equals(BigInteger.ONE)) 
				return false;
			if(x.equals(n_minus_one))
				return true;
		}
		return false;
    }
	//Square-and-Multiply algorithm to calculate (x^e) mod m
	static BigInteger sqm(BigInteger x, BigInteger e, BigInteger m){
		BigInteger r = new BigInteger("1");
		while(e.compareTo(BigInteger.ZERO)>0)
		{
			if(e.mod(big_two).equals(BigInteger.ONE))//if e is odd
				r = r.multiply(x).mod(m);//(r*x) mod m
			e = e.divide(big_two);
			x = x.multiply(x).mod(m);
		}	
		return r;
	}
}