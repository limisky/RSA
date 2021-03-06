import java.io.*;
import java.math.BigInteger;
import java.util.Hashtable;
import java.util.Random;

public class RSA {
	static int L = 2;					//size of block - the encryption deals with L characters at a time
	final static int PQ_BITLENGTH = 32;		//size of p,q , the maximum allowed number is 256 for this lab, minimum is 16
	final static int KEY_BITLENGTH = 32;	//size of e, can not be greater then 2*PQ_BITLENGTH
	
	
	// Have to increase heap space
	// http://edwards.sdsu.edu/research/index.php/daniel/236-increasing-heap-size-in-eclipse
	
	static Random rnd = new Random();
	static BigInteger big_two = new BigInteger("2"); //BigInteger is Java's way of using big integer. This declares two.
	
	public static void main(String[] args) {
		int r;	
		String strLine;
		int groupNr;
		double totalTime = 0;
		long endTime;
		long startTime;
		double time;

		/* This will solve the groups crypto files*/ 
  		try{
  			File statistic = new File("lab2_statstics.txt");
  			for (int i = 8; i <= 8; i++) {
  				if (i == 5 || i == 7)
  					continue;
				L = 1;
  				groupNr = i;
				for (int y = 2; y <= 2; y++) {
					L = y;
					r = L*8 -1;
					for (int j = 2; j <= 2; j++) {
						// r++;
						FileWriter statisticWriter = new FileWriter(statistic, true);
						System.out.println("STARTING BREAKING CIPHER \n");
						System.out.println("-----------------------------\n");
						System.out.println("Group: "+groupNr+"\n");
						System.out.println("L: "+L+"\n");
						System.out.println("r: "+r+"\n");
 	 		
						FileInputStream fstream = new FileInputStream("rsa_group"+groupNr+"_"+L+".pub");
  						DataInputStream in = new DataInputStream(fstream);
  						BufferedReader br = new BufferedReader(new InputStreamReader(in));

  						strLine = br.readLine();
  						BigInteger publicKey_E = new BigInteger(strLine);
  						strLine = br.readLine();
						BigInteger publicKey_N = new BigInteger(strLine);
						in.close();

  						fstream = new FileInputStream("rsa_group"+groupNr+"_"+L+".crypto");
  						in = new DataInputStream(fstream);
  						br = new BufferedReader(new InputStreamReader(in));
  						File crypto = new File("rsa_group"+groupNr+"_"+L+"_"+r+".plain");
  					
						crypto.createNewFile();
										
						startTime = System.nanoTime();

  						while ((strLine = br.readLine()) != null)   {
  							FileWriter cryptoWriter = new FileWriter(crypto, true);

  							BigInteger cipher = new BigInteger(strLine);
  							Object[] hashtables = generatePossibleCiphers(publicKey_E,  publicKey_N, r);
							@SuppressWarnings("unchecked")
							Hashtable<BigInteger,BigInteger> hashtableIndex = (Hashtable<BigInteger,BigInteger>) hashtables[1];
							@SuppressWarnings("unchecked")
							Hashtable<BigInteger,BigInteger> hashtablePowMod = (Hashtable<BigInteger,BigInteger>) hashtables[0];
		
							Hashtable<BigInteger,BigInteger> hashtableInverses = generateInverses(hashtableIndex, r, publicKey_N);

							BigInteger plain = breakCipher(hashtablePowMod, hashtableInverses, publicKey_N, cipher, r);
							if (L == 1){
								cryptoWriter.write((char) plain.intValue());
							}
							else if (L == 2) {
								int[] plain_temp = bitShift(plain, L);
								cryptoWriter.write((char) plain_temp[0]);
								cryptoWriter.write((char) plain_temp[1]);
							}
							else if (L == 3) {
								int[] plain_temp = bitShift(plain, L);
								cryptoWriter.write((char) plain_temp[0]);
								cryptoWriter.write((char) plain_temp[1]);
								cryptoWriter.write((char) plain_temp[2]);
							}
							
							cryptoWriter.close();
  						}

  						endTime = System.nanoTime();
  						in.close();
  					
  						time = ((endTime - startTime) / 1000000000);
  						totalTime += (time / 60)/60;
  						statistic.createNewFile();
  						statisticWriter.write(groupNr + " & " + r + " & " + L + " & " + time + "\n");
  						statisticWriter.close();
  					
  						System.out.println("-----------------------------\n");
 	 					System.out.println("FINISHED BREAKING CIPHER\n");
 	 					System.out.println("Total Time: "+totalTime+"\n");
					}
				}
  			}
  			
    	} catch (Exception e){ 
  			System.err.println("Error: " + e.getMessage());
 	 	}
 	 	/* This will solve the assistance crypto files */ 
 	 	try{
  			File statistic = new File("lab2_assist_statstics.txt");
  			

  			for (int i = 1; i <= 1; i++) {
  				if (i == 5 || i == 7)
  					continue;
				L = 1;
  				groupNr = i;
				for (int y = 3; y <= 7; y++) {
					L = y;
					r = L*8;
					FileWriter statisticWriter = new FileWriter(statistic, true);
					System.out.println("STARTING BREAKING ASSISTANT CIPHER \n");
					System.out.println("-------------------------------------\n");
					System.out.println("L: "+L+"\n");
					System.out.println("r: "+r+"\n");
 	 		
					FileInputStream fstream = new FileInputStream("rsa_assist_"+L+".pub");
  					DataInputStream in = new DataInputStream(fstream);
  					BufferedReader br = new BufferedReader(new InputStreamReader(in));

  					strLine = br.readLine();
  					BigInteger publicKey_E = new BigInteger(strLine);
  					strLine = br.readLine();
					BigInteger publicKey_N = new BigInteger(strLine);
					in.close();

  					fstream = new FileInputStream("rsa_assist_"+L+".crypto");
  					in = new DataInputStream(fstream);
  					br = new BufferedReader(new InputStreamReader(in));
  					File crypto = new File("rsa_assist_"+L+".plain");
  			
					crypto.createNewFile();
					
					startTime = System.nanoTime();

  					while ((strLine = br.readLine()) != null)   {
  						FileWriter cryptoWriter = new FileWriter(crypto, true);
  						BigInteger cipher = new BigInteger(strLine);
  						Object[] hashtables = generatePossibleCiphers(publicKey_E,  publicKey_N, r);
						@SuppressWarnings("unchecked")
						Hashtable<BigInteger,BigInteger> hashtableIndex = (Hashtable<BigInteger,BigInteger>) hashtables[1];
						@SuppressWarnings("unchecked")
						Hashtable<BigInteger,BigInteger> hashtablePowMod = (Hashtable<BigInteger,BigInteger>) hashtables[0];
		
						Hashtable<BigInteger,BigInteger> hashtableInverses = generateInverses(hashtableIndex, r, publicKey_N);

						BigInteger plain = breakCipher(hashtablePowMod, hashtableInverses, publicKey_N, cipher, r);
						if (L == 1){
								cryptoWriter.write((char) plain.intValue());
							}
							else if (L == 2) {
								int[] plain_temp = bitShift(plain, L);
								cryptoWriter.write((char) plain_temp[0]);
								cryptoWriter.write((char) plain_temp[1]);
							}
							else if (L == 3) {
								int[] plain_temp = bitShift(plain, L);
								cryptoWriter.write((char) plain_temp[0]);
								cryptoWriter.write((char) plain_temp[1]);
								cryptoWriter.write((char) plain_temp[2]);
							}
							
							cryptoWriter.close();
  						}

  					endTime = System.nanoTime();
  					in.close();
  					
  					time = ((endTime - startTime) / 1000000000);
  					totalTime += (time/60)/60;
  					statistic.createNewFile();
  					statisticWriter.write(groupNr + " & " + r + " & " + L + " & " + time + "\n");
  					statisticWriter.close();
  					
  					System.out.println("-----------------------------\n");
 	 				System.out.println("FINISHED BREAKING CIPHER\n");
 	 				System.out.println("Total Time: "+totalTime+"\n");

				}
  			}
  			

  			
    	} catch (Exception e){ 
  			System.err.println("Error: " + e.getMessage());
 	 	}
	}
	
	static void createRSACipher(){
				
		//start time capture
		long startTime = System.nanoTime();
		
		//Generates the 3 keys that RSA need.
		//KeyPair[0] = e
		//KeyPair[1] = n
		//KeyPair[2] = d
		//KeyPair[3] = n
		BigInteger[] keyPair = generateKey();
		
		
		//reads our file
		int[] plainText = readFile("rsa_group9_"+L+".plain");
		
		//Encrypt plaintext
		BigInteger[] cipher = encrypt(plainText,L,keyPair[0],keyPair[1]);

		//end of time capture
		long endTime = System.nanoTime();
		System.out.println("Run time: "+(endTime-startTime)+" ns"); 
		
		//Writes files :
		// "rsa_group9_"+L+".crypto"
		// "rsa_group9_"+L+".key"
		// "rsa_group9_"+L+".pub"
		writeFile(cipher,keyPair);
	}
	
	
	//c * inv(hashtableIndex[i], n) mod n = X
			//if j = hashtablePowMod[X] exist, then
			// m = i * j mod n
	
		// or -1 if not found
	static BigInteger breakCipher(Hashtable<BigInteger, BigInteger> hashtablePowMod,
			Hashtable<BigInteger, BigInteger> hashtableInverses, BigInteger publicKey_N, BigInteger cipher, int r)
	{
		BigInteger Max = big_two.pow(r);
		
		//is the multiplication before or after modulo?
		for (BigInteger i = Max; i.compareTo(BigInteger.ZERO) == 1 ; i = i.subtract(BigInteger.ONE)){
			
			BigInteger X = cipher.multiply(hashtableInverses.get(i).mod(publicKey_N));
			BigInteger j = hashtablePowMod.get(X);
			BigInteger minusOne = new BigInteger("-1");

				//if X was in the table then we found our plaintext
				if (j != null){
					//Found m!!
					return  i.multiply(j).mod(publicKey_N);
				}
		}
		return new BigInteger("-1");
	}
	
	//inv(i^e mod n, n)
	static Hashtable<BigInteger, BigInteger> generateInverses(Hashtable<BigInteger,BigInteger> hashtableIndex,int r , BigInteger publicKey_N)
	{
		//Make a hashtable of inv(hashtableIndex[i], n) with index i from 0 to 2^r 
		Hashtable<BigInteger, BigInteger> hashtableInverses = new Hashtable<BigInteger, BigInteger>();
		BigInteger Max = big_two.pow(r);
		
		for (BigInteger i = Max; i.compareTo(BigInteger.ZERO) == 1 ; i = i.subtract(BigInteger.ONE)){
			BigInteger Temp = inv(hashtableIndex.get(i), publicKey_N);
			hashtableInverses.put(i, Temp);		
		}
		return hashtableInverses;
	}
	
	
	//for i := 1 to 2^r
		//	table[i] = i^e mod n
	
	static Object[] generatePossibleCiphers(BigInteger publicKey_E, BigInteger publicKey_N,int r)
	{
		Hashtable<BigInteger, BigInteger> hashtablePowMod = new Hashtable<BigInteger, BigInteger>();
		Hashtable<BigInteger, BigInteger> hashtableIndex = new Hashtable<BigInteger, BigInteger>();	
		
		BigInteger Max = big_two.pow(r);
		BigInteger Temp;
		
		// Do the loop 2^r times. Beginning from 2^r and ending on 0
		for (BigInteger i = Max; i.compareTo(BigInteger.ZERO) == 1 ; i = i.subtract(BigInteger.ONE)){
			Temp = sqm(i,publicKey_E,publicKey_N);
			hashtablePowMod.put(Temp, i);
			hashtableIndex.put(i, Temp);
			//System.out.println(Temp);
		}
		
		//Hashtable<BigInteger, BigInteger> hashtableReturn[] = new Hashtable<BigInteger, BigInteger>[2];
		//return hashtableIndex;
		return new Object[]{hashtablePowMod, hashtableIndex};
	}
	
	static int[] bitShift(BigInteger crypto, int L) {
		int[] plain = new int[L];
		BigInteger mask= new BigInteger("255"); //used for masking out bits

		for (int j = L-1;j>=0;j--)
			{
				if (L==3){
					if (j==0)			
						plain[(L-j-1)] = crypto.and(mask).intValue();
					else if(j==1)
						plain[(L-j-1)] = crypto.shiftRight(8).and(mask).intValue();
					else //(j==2)
						plain[(L-j-1)] = crypto.shiftRight(16).intValue();
				}
				else if (L==2){
					if (j==0)
						plain[(L-j-1)] = crypto.and(mask).intValue();
					else //(j==1)
						plain[(L-j-1)] = crypto.shiftRight(8).intValue();
				}
				else // L==1
					plain[(L-j-1)] = crypto.intValue();
		}
		return plain;

	}

	
	

	
	/**
	 * Writes cipher, private and public key to files.
	 * 
	 * Write cipher into file "rsa_group9_L.crypto"
	 * Write private key into file "rsa_group9_L.key"
	 * Write public key into file "rsa_group9_L.pub"
	 * 
	 * @param cipher  - The cipher
	 * @param keyPair - Public and private keys
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
	
	/**
	 * d(k) = x^d mod n where x represent L characters
	 * 
	 * @param cipher - cipher
	 * @param L 	 - Chars to be decrypted at the same time.
	 * @param key_d  - d of public key pair in RSA
	 * @param key_n  - n of public key pair in RSA
	 * @return plaintext, each char on a own space
	 */
	static int[] decrypt(BigInteger[] cipher, int L, BigInteger key_d, BigInteger key_n)
	{
		int[] plain = new int[102];
		BigInteger tempCipher;
		BigInteger tempPlain;
		BigInteger mask= new BigInteger("255"); //used for masking out bits
		for(int i=0;cipher[i]!=null;i++)
		{
			tempCipher = cipher[i];
			
			tempPlain = sqm(new BigInteger(""+ tempCipher),key_d,key_n);//decrypt		
			
			for (int j = L-1;j>=0;j--)
			{
				if (L==3){
					if (j==0)			
						plain[L*i+(L-j-1)] = tempPlain.and(mask).intValue();
					else if(j==1)
						plain[L*i+(L-j-1)] = tempPlain.shiftRight(8).and(mask).intValue();
					else //(j==2)
						plain[L*i+(L-j-1)] = tempPlain.shiftRight(16).intValue();
				}
				else if (L==2){
					if (j==0)
						plain[L*i+(L-j-1)] = tempPlain.and(mask).intValue();
					else //(j==1)
						plain[L*i+(L-j-1)] = tempPlain.shiftRight(8).intValue();
				}
				else // L==1
					plain[L*i+(L-j-1)] = tempPlain.intValue();
					
			
				
			}
		}
		return plain;
	}
	
	/**
	 * e(k) = x^e mod n where x represent L characters
	 * 
	 * @param plain - plaintext
	 * @param L 	 - Chars to be encrypted at the same time.
	 * @param key_e  - e of public key pair in RSA
	 * @param key_n  - n of public key pair in RSA
	 * @return encrypted plaintext with RSA
	 */
	static BigInteger[] encrypt(int[] plain, int L, BigInteger key_e, BigInteger key_n)
	{
		BigInteger[] cipher= new BigInteger[100/L+1];
		int temp = 0;
		int index = 0;
		//end of plaintext is -1
		for(int i=0 ,l=L ;plain[i]!=-1 ;i++)
		{
			//encode L characters at a time in 8 bits
			temp = (temp << 8) + plain[i];
			if(l==1)
			{
				cipher[index] = sqm(new BigInteger(temp+""),key_e,key_n);//encrypt
				l=L;
				index++;
				temp = 0;
			}
			else
				l--;
		}
		return cipher;
	}
	
	/**
	 * Read plain text into an integer array,
	 * characters are encoded on 8 bits (typically Latin1 encoding)
	 * 
	 * @param fileName - file to read from
	 * @return the first 102 chars, one char in one space
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
		//----------------We pad the message with zeroes-------------
		//------------------------------------------------------------
		while(i%L!=0)
		{
			plain[i]=0;
			i++;
		}
        return plain;
	}
	/**
	 * Generate key pair for RSA
	 * @return 
	 * 		result[0]&result[1] - public key (e,n)
	 * 		result[2]&result[3] - private key (d,n)
	 */
	static BigInteger[] generateKey()
	{
		BigInteger result[] = new BigInteger[4];
		BigInteger p = generatePrime(PQ_BITLENGTH);
		BigInteger q,phi,n,e,d;
		int count = 0;
		do{
			q = generatePrime(PQ_BITLENGTH);
		}while (q.equals(p));
		n = p.multiply(q);//n = p*q
		
		phi = (p.subtract(BigInteger.ONE)).multiply((q.subtract(BigInteger.ONE)));//phi = (p-1)*(q-1)
		do{
			e = new BigInteger(KEY_BITLENGTH,rnd);
			e = e.setBit(KEY_BITLENGTH-1); //make sure e uses KEY_BITLENGTH bits
			count++;
		}while(!gcd(e,phi).equals(BigInteger.ONE)||e.compareTo(phi)>=0||count<10000);
		//Try a maximumof 10 000 times before using new probable prime numbers
		if(count == 10000)
		{
			return generateKey();
		}		

		d = inv(e,phi);
		
		result[0]=e;
		result[1]=n;
		result[2]=d;
		result[3]=n;
		
		return result;
	}

	/**
	 * Use Extended Euclidean Algorithm to find y = inv(x,n) is such that (x*y=1) mod n
	 * ax+by=1 return x+b
	 * 
	 * @param x - 
	 * @param n - 
	 * 
	 * @return x+b
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
	
	
	/**
	 * calculate the greatest common divisor of a and b
	 * 
	 * @param a - 
	 * @param b - 
	 * 
	 * @return gcd(a,b)
	 */
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
	
	/**
	 * generate a strong probable prime with the specified bitLength
	 * 
	 * @param bitLength - number of bits the prime should be
	 * @return a probable prime
	 */
	static BigInteger generatePrime(int bitLength){
		BigInteger x = new BigInteger(bitLength,rnd);//[0,2^bitLength-1]
		//change the most significant bit to 1 to make sure it meets the bitLength requirement
		x = x.setBit(bitLength-1);
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
    
    /**
	 * //Square-and-Multiply algorithm to calculate (x^e) mod m with BigIntegers
	 * 
	 * @param x - 
	 * @param e - 
	 * @param m - 
	 * @return (x^e) mod m
	 */
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