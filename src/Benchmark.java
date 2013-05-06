import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.math.BigInteger;

public class Benchmark {
	public static void main(String[] args) {
		int p, q, e, d, L, text_length;
		int[] plaintext;
		File encrypt = new File("benchmark_encrypt");
		File decrypt = new File("benchmark_decrypt");
		File generateKey = new File("benchmark_generatekey");
		BigInteger[] keyPair;
		BigInteger[] cipher;
		long startTime;
		long endTime;
		long time;
		try {
			encrypt.createNewFile();
			decrypt.createNewFile();
			generateKey.createNewFile();
		} catch(Exception a) {
			a.printStackTrace();
		}
		
		// benchmark for generateKey()
		startTime = System.nanoTime();
		keyPair = generateKey();
		endTime = System.nanoTime();
		time = endTime - startTime;
		writeFile(generateKey, 1, 2, 3, 4, 5, 6, 100);

		// bechmark for encrypt()
		startTime = System.nanoTime();
		cipher = encrypt(plaintext, L, keyPair[0], keyPair[1]);
		endTime = System.nanoTime();
		time = endTime - startTime;
		writeFile(encrypt, 1, 2, 3, 4, 5, 6, 100);

		// benchmark for decrypt()
		startTime = System.nanoTime();
		keyPair = decrypt();
		endTime = System.nanoTime();
		time = endTime - startTime;
		writeFile(decrypt, 1, 2, 3, 4, 5, 6, 100);
	}

	static void
	writeStatistics(FileWriter file, File fileName, int p, int q,
		int e, int d, int L, int text_length, int time) {
		try {
			file.write(
			p + " & " + 
			q + " & " +
			e + " & " +
			d + " & " +
			L + " & " +
			text_length + " & " +
			time + "\n");
		} catch(Exception a) {
			a.printStackTrace();
		}	
	}

	static void
	writeFile(File fileName, int p, int q, int e, int d,
		int L, int text_length, int time) {
		try {
			FileWriter statistic = new FileWriter(fileName);
			for (int i = 0; i < 10; i++) 
				writeStatistics(statistic, fileName, 1, 2, 3, 4, 5, 6, 100);
			statistic.close();
		} catch(Exception a) {
			a.printStackTrace();
		}
	}

}












