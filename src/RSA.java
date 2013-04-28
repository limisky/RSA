import java.math.BigInteger;
import java.util.Random;

public class RSA {
	final static int PRIME_BITLENGTH = 16;
	final static BigInteger E = new BigInteger("65537");
	static Random rnd = new Random();
	static BigInteger big_two = new BigInteger("2");
	
	public static void main(String[] args) {
		
		System.out.println(sqm(new BigInteger("3"),new BigInteger("5"),new BigInteger("7")).toString());
		
		BigInteger testprime = new BigInteger("11111111111111111111111");
		System.out.println(miller_rabin(testprime,1) ? "PRIME" : "COMPOSITE");
		
		System.out.println(generatePrime(PRIME_BITLENGTH));
		
		BigInteger[] keyPair = generateKey();
		System.out.println("public key: (" + keyPair[0] + "," + keyPair[1] + ")");
		System.out.println("private key: (" + keyPair[2] + "," + keyPair[3] + ")");
		
		
		
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
		BigInteger p = generatePrime(PRIME_BITLENGTH);
		BigInteger q;
		BigInteger phi;
		BigInteger n;
		BigInteger d;
		do{
			do{
				q = generatePrime(PRIME_BITLENGTH);
			}while (q.equals(p));
			n = p.multiply(q);
			//phi = (p-1)*(q-1)
			phi = (p.subtract(BigInteger.ONE)).multiply((q.subtract(BigInteger.ONE)));
		}while(!gcd(E,phi).equals(BigInteger.ONE)||E.compareTo(phi)>=0);
		d=inv(E,phi);
		
		result[0]=E;
		result[1]=n;
		result[2]=d;
		result[3]=n;
		
		return result;
	}
	//Use Extended Euclidean Algorithm to find y = inv(x,n) is such that (x*y=1) mod n
	static BigInteger inv(BigInteger x, BigInteger n)
	{
		/*
		 * Something here
		 */
		
		return null;
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
		x = x.or(BigInteger.ONE.shiftLeft(bitLength-1));//[2^(bitLength-1),2^bitLength-1]
		if(x.mod(big_two).equals(BigInteger.ZERO))//if n is even
			x = x.add(BigInteger.ONE);
		while(!miller_rabin(x,1))
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

