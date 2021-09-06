package count.machine;

public class Connection
{
	public Connection(Machine upstream, Machine downstream, int width)
	{
		this.upstream = upstream;
		this.downstream = downstream;
		this.width = width;
	}
	
	Machine upstream;
	Machine downstream;
	int width;
	

}
